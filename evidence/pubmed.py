"""Read-only NCBI PubMed retrieval. Bibliographic verification is not method validation."""
import json, re, threading, time, urllib.error, urllib.parse, urllib.request
import xml.etree.ElementTree as ET

_lock = threading.Lock()
_last = 0.0
BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

def get(endpoint, params):
    global _last
    params = dict(params, tool='KilnBioEvidence')
    url = BASE + endpoint + '?' + urllib.parse.urlencode(params)
    for attempt in range(3):
        with _lock:
            time.sleep(max(0, .4 - (time.monotonic() - _last)))
            _last = time.monotonic()
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'KilnBio/0.2 evidence audit'})
            with urllib.request.urlopen(req, timeout=40) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError):
            if attempt == 2: raise
            time.sleep(2 ** attempt)

def search(query, limit=10, offset=0):
    if not isinstance(query, str) or not query.strip() or len(query) > 3000:
        raise ValueError('Provide a nonempty public-literature query of at most 3000 characters')
    if type(limit) is not int or not 1 <= limit <= 50:
        raise ValueError('limit must be 1..50')
    if type(offset) is not int or not 0 <= offset <= 9950:
        raise ValueError('offset must be 0..9950; partition queries exceeding PubMed search limits')
    data = json.loads(get('esearch.fcgi', {'db':'pubmed', 'term':query, 'retmode':'json', 'retmax':limit, 'retstart':offset}))
    result = data['esearchresult']
    return {'query':query, 'total':int(result['count']), 'offset':offset, 'pmids':result['idlist'],
            'query_translation':result.get('querytranslation', '')}

def content(node):
    return ''.join(node.itertext()).strip() if node is not None else ''

def parse_records(xml, include_abstract=False):
    root = ET.fromstring(xml)
    records = []
    for article in root.findall('PubmedArticle'):
        citation = article.find('MedlineCitation')
        body = citation.find('Article')
        ids = {n.attrib.get('IdType'):content(n) for n in article.findall('./PubmedData/ArticleIdList/ArticleId')}
        links = [{'type':n.attrib.get('RefType'), 'pmid':content(n.find('PMID')),
                  'citation':content(n.find('RefSource'))}
                 for n in citation.findall('./CommentsCorrectionsList/CommentsCorrections')]
        types = [content(n) for n in body.findall('./PublicationTypeList/PublicationType')]
        flags = [p for p in types if p in ('Retracted Publication','Retraction of Publication','Expression of Concern','Preprint')]
        flags += [n['type'] for n in links if n['type'] in ('RetractionIn','ExpressionOfConcernIn')]
        record = {'pmid':content(citation.find('PMID')), 'doi':ids.get('doi'), 'pmcid':ids.get('pmc'),
                  'title':content(body.find('ArticleTitle')), 'journal':content(body.find('./Journal/Title')),
                  'year':content(body.find('./Journal/JournalIssue/PubDate/Year')),
                  'article_dates':[{'type':d.attrib.get('DateType'), 'year':content(d.find('Year')),
                                    'month':content(d.find('Month')), 'day':content(d.find('Day'))}
                                   for d in body.findall('ArticleDate')],
                  'publication_types':types, 'related_notices':links, 'safety_flags':sorted(set(flags)),
                  'authors':[content(a.find('LastName')) + ' ' + content(a.find('Initials'))
                             for a in body.findall('./AuthorList/Author')],
                  'url':'https://pubmed.ncbi.nlm.nih.gov/' + content(citation.find('PMID')) + '/'}
        if include_abstract:
            record['abstract'] = '\n'.join(content(n) for n in body.findall('./Abstract/AbstractText'))
        records.append(record)
    return records

def fetch(pmids, include_abstract=False):
    if not isinstance(pmids, list) or not 1 <= len(pmids) <= 50 or any(not re.fullmatch(r'[1-9][0-9]{0,9}', str(p)) for p in pmids):
        raise ValueError('Provide 1..50 numeric PubMed identifiers')
    return parse_records(get('efetch.fcgi', {'db':'pubmed','id':','.join(map(str,pmids)), 'retmode':'xml'}), include_abstract)
