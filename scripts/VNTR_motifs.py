#!/usr/bin/env python

from collections import Counter
import regex
import sys


TR = sys.argv[1]
error_rate = int(0.2 * len(TR))
print(regex.findall(f"(?e)({TR}){{e<4}}", sys.argv[2]))
print(Counter(regex.findall(f"(?e)({TR}){{e<{error_rate}}}", sys.argv[2])).most_common(5))
