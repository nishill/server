import sys
import os
import re
def main():
  rep = {"source_name": "sourceName", "source_version": "sourceVersion"}
  rep = dict((re.escape(k), v) for k, v in rep.iteritems())
  pattern = re.compile("|".join(rep.keys()))
  with open ("new_out.json", "wt" ) as fout:
    with open ("out.json", "rt" ) as fin:
      for line in fin:
        fout.write(pattern.sub(lambda m: rep[re.escape(m.group(0))], line))
		

if __name__ == "__main__":
    main()
