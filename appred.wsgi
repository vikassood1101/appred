#! /home/ubuntu/appred/env/bin python3.8
import site
site.addsitedir('/var/www/html/appred/env/lib/python3.8/site-packages')
import sys
sys.path.insert(0, '/var/www/html/appred')

from appred import app as application
application.secret_key = 'gkuotkgjbnbbfpptiuxvznametyworghghvdkkdfllorrowgcbry'