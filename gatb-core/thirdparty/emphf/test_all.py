#!/usr/bin/env python

import os
import sys
from subprocess import check_call

DEFAULT_FILE = '/usr/share/dict/words'

EXES = [
    ('compute_mphf_seq', 'test_mphf'),
    ('compute_mphf_scan', 'test_mphf'),
    ('compute_mphf_scan_mmap', 'test_mphf'),
    ('compute_mphf_hem', 'test_mphf_hem'),
    ]

def main(argv):
    if len(argv) == 1:
        filename = DEFAULT_FILE
        print >> sys.stderr, "Using default file %s" % filename
        print >> sys.stderr, "To use another file:"
        print >> sys.stderr, "\t%s <filename>" % argv[0]
    else:
        filename = argv[1]
        print >> sys.stderr, "Using default file %s" % filename


    for constructor, tester in EXES:
        print >> sys.stderr
        print >> sys.stderr, '=' * 4, 'Testing %s' % constructor, '=' * 40
        mphf_name = 'mphf.output.bin'

        check_call(['./' + constructor, filename, mphf_name])
        check_call(['./' + tester, filename, mphf_name, '--check'])
        check_call(['rm', mphf_name])


if __name__ == '__main__':
    main(sys.argv)
