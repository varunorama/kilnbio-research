"""Static regression guard: no general-purpose LLM integrations in shipped code.

This catches known imports, endpoints and CLI launch patterns. It is not a
network sandbox or a proof about every transitive third-party dependency.
"""
import ast
from pathlib import Path
import re
import sys

FORBIDDEN_IMPORTS = {'openai','anthropic','ollama','litellm','llama_cpp',
                     'langchain','langchain_openai','langchain_anthropic',
                     'google.generativeai','google.genai'}
ENDPOINTS = re.compile(r'api\.openai\.com|api\.anthropic\.com|generativelanguage\.googleapis\.com|'
                       r'bedrock-runtime\.[a-z0-9-]+\.amazonaws\.com', re.I)
COMMANDS = re.compile(r'(?:\b(?:exec|command|subprocess\.(?:run|Popen|call))\b[^\n]{0,160})'
                      r'[\s\x22\x27](?:codex|claude|ollama)(?:[\s\x22\x27]|$)',re.I)


def inspect(path):
    text=path.read_text(encoding='utf-8-sig')
    problems=[]
    if ENDPOINTS.search(text): problems.append('general-purpose model endpoint')
    if COMMANDS.search(text): problems.append('general-purpose model CLI')
    if path.suffix=='.py':
        for node in ast.walk(ast.parse(text,filename=str(path))):
            names=([n.name for n in node.names] if isinstance(node,ast.Import)
                   else [node.module or ''] if isinstance(node,ast.ImportFrom) else [])
            for name in names:
                if any(name==bad or name.startswith(bad+'.') for bad in FORBIDDEN_IMPORTS):
                    problems.append('general-purpose model import: '+name)
    return problems


def scan(root):
    results=[]
    for path in root.rglob('*'):
        rel=path.relative_to(root)
        if any(x.startswith('.') or x in {'tests','ci','__pycache__'} for x in rel.parts): continue
        if path.suffix.lower() not in {'.py','.r','.sh','.nf'}: continue
        for reason in inspect(path): results.append((str(rel),reason))
    return results


if __name__=='__main__':
    found=scan(Path(__file__).resolve().parents[1])
    for path,reason in found: print(path+': '+reason)
    print('Runtime source policy: '+('FAIL' if found else 'PASS (static checks only)'))
    sys.exit(bool(found))
