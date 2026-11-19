"""Guard GitHub math delimiters and TeX grouping in repository documentation."""
from pathlib import Path
import re
import unittest

ROOT=Path(__file__).resolve().parents[1]


class MarkdownMath(unittest.TestCase):
    def test_math_blocks_and_inline_groups(self):
        for path in ROOT.rglob('*.md'):
            if '.git' in path.parts:continue
            text=path.read_text(encoding='utf-8')
            blocks=re.findall(r'(?m)^```math\s*\n(.*?)\n```',text,re.S)
            self.assertEqual(len(blocks),len(re.findall(r'(?m)^```math\s*$',text)),str(path))
            prose=re.sub(r'(?m)^```[^\n]*\n.*?^```[^\n]*(?:\n|$)','',text,flags=re.S)
            self.assertNotRegex(prose,r'(?m)^\$\$\s*$',str(path))
            expressions=blocks+re.findall(r'\$`([^`]+)`\$',prose)
            for tex in expressions:
                self.assertNotIn('\\operatorname',tex,str(path))
                self.assertNotIn('<',tex,'Use \\lt to avoid GitHub HTML parsing: '+str(path))
                depth=0
                for brace in re.findall(r'(?<!\\)[{}]',tex):
                    depth+=1 if brace=='{' else -1
                    self.assertGreaterEqual(depth,0,str(path)+' '+tex)
                self.assertEqual(depth,0,str(path)+' '+tex)


if __name__=='__main__':unittest.main()
