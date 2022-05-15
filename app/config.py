from distutils.debug import DEBUG
import os
basedir = os.path.abspath(os.path.dirname(__file__))

class Config():
    DEBUG = False
    
class LocalDevelopementConfig(Config):
    DEBUG = True
    SECRET_KEY = "324hretrgbrr4hh3jhhvh43jre" # Strong, random, difficult key

class ProductionDevelopementConfig(Config):
    DEBUG = False
    SECRET_KEY = "u23brn234sdfasddfmadsandalv23sdf" # Strong, random, difficult key