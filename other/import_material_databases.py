
import os,sys
import urllib, urllib2, urlparse

from httplib import HTTPConnection
from sgmllib import SGMLParser

#Custom class to just locate url links
class URLMaterialFinder(SGMLParser):
    def load_url(self, urladdress):
        conn = urllib2.urlopen(urladdress)
        self.feed(conn.read())
        conn.close()

    def reset(self):
        SGMLParser.reset(self)
        self.urls = []

    def start_a(self, attrs):
        for k, v in attrs:
            if k=='href' and v.startswith("Material"):
                self.urls.append(v)

    def download_urls(self, baseurl, dirname):
        try:
            os.mkdir(dirname)
        except:
            pass
        
        #Download all urls in urllist        
        if os.path.exists(dirname) and os.path.isdir(dirname):
            for url in self.urls:
                _,filename = os.path.split(url)
                print "Downloading file", filename
                try:
                    urllib.urlretrieve(urlparse.urljoin(baseurl,url), os.path.join(dirname,filename))
                except:
                    print "Error: Couldn't download file!"

#Request data from luxpop
print "Downloading list of LUXPOP materials"
materialurls =  URLMaterialFinder()
materialurls.load_url("http://luxpop.com/RefractiveIndexList.html")
materialurls.close()

#Download material files
print "Downloading LUXPOP material files"
materialurls.download_urls("http://luxpop.com/",os.path.join(os.path.pardir, "Polymode/data/luxpop"))

print "The data from SOPRA can also be downloaded at the location:"
print "http://www.sopra-sa.com/index2.php?goto=dl&rub=4"

