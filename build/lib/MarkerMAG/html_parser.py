from requests import get
from bs4 import BeautifulSoup
from html.parser import HTMLParser

html_file = '/Users/songweizhi/Desktop/vis_msa/combined_MView_tmp.html'

soup = BeautifulSoup(open(html_file), features='html.parser')


print(soup.body)


print(soup.body.get_text())

