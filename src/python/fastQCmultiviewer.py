#! /usr/bin/python

def keyCheck(key, dic):
	try:
		dic[key]
		return True
	except KeyError:
		return False

def createIframe(html):
	string = "<td><iframe src='"+html+"' width='1500px' height='10000px' scrolling='no'></iframe></td>"
	return string

def siteHtml(images):
	html = """
			<!DOCTYPE html>
				<head></head>

				<body>
				<table>
				<tr>
				"""+"".join(images)+"""
				</tr>
				</table>

				</body>

				</html>
	"""
	return html

import sys
import os

inputs = sys.argv

inputs.pop(0)
inDic = {}
for inp in range(len(inputs)-1):
	if inputs[inp][0] == "-":
		inDic[inputs[inp][1:]] = inputs[inp+1]

outFile = open(inDic['out_file'], "w")
limit = [0,10]
if keyCheck("range",inDic):
	limit = inDic["range"].strip().split(",")


table = []
if keyCheck("samples", inDic):
	samples = inDic["samples"].strip.split(",")
if keyCheck("in_dir", inDic):
	rows = 0
	folders = os.walk(inDic["in_dir"])
	for root, subFolders, files in folders:
		#print root
		folders = root.split("/")
		relative_folder = folders[-2]+"/"+folders[-1]
		for f in files:
			if f[-4:] == "html":
				table.append(createIframe(relative_folder+"/"+f))
				if rows > int(limit[1]):
					break
				rows += 1
	table.sort()
	outFile.write(siteHtml(table))