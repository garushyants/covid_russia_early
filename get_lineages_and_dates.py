#!/usr/bin/python3

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace
import random
import re
from datetime import datetime

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='', type='str')
#parser.add_option('-l', '--leafnames', default= '',help='', type='str')
parser.add_option('-t', '--nonbinary', default= 1, help='', type='int')
#parser.add_option('-r', '--root', default= "hCoV-19/Wuhan/WH01/2019|EPI_ISL_406798|2019-12-26", help='', type='str')
parser.add_option('-c','--cluster', default= '', help='',type='str')
parser.add_option('-b','--bootstrap',default = 1, type = 'int')
parser.add_option('-v','--verbose',help='print additional info on the screen',action='store_true')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.bootstrap) #input tree is from treetime and nonbinary
#root = options.root
nbchecker = options.nonbinary
#clid=options.number
#if options.leafnames == '':
#	print("No file with leaves names provided")	
###
cutoff = 1e-05

####
#reroot
#tree.set_outgroup(tree&root) #the rerooting parameter was removed, because Ksusha already rerooted the trees
####
#Read correct names
leaf_dates = dict()
leaf_names = dict()

##
def rename_leaves(node):
	for ll in node:
		l = ll.name
		#print(l)
		pattern = re.compile("hCoV-19/", re.IGNORECASE)
		prename = pattern.sub("", l)
		pnum = prename.count('|')
		name = ""
		gisid = ""
		datestr = ""
		num = 0
		if pnum ==2:
			name,gisid,datestr = prename.split("|")
		if pnum ==3:
			name,gisid,datestr,num = prename.split("|")
		datestr_corr = datestr
		if ("?" in datestr) or (not datestr.startswith("20")):
			datestr_corr = "2300-01-01"
		datecheck = datestr_corr.split("-")
		if len(datecheck)<3:
			if len(datecheck) ==2:
				datestr_corr = datestr_corr+"-28"#set to max possible date
			if len(datecheck) ==1:
				datestr_corr = datestr_corr+"-12-31"#set to the latest possible date if data is missing
		#print(datestr_corr)
		date = datetime.strptime(datestr_corr, "%Y-%m-%d")
		newdateformat = date.strftime("%b-%d")
		newname=""
		if "Russia" in l or "FMBA" in l:
			namedate = str(newdateformat)
			if ("?" in datestr) or (not datestr.startswith("20")):
				namedate = "Date unknown"
			elif (len(datecheck) < 3): 
				namedate = datestr
			#print(datestr)
			#print(namedate)
			newname = name+"|"+namedate
		else:
			newname = name+"|"+gisid+"|"+str(newdateformat)
		leaf_names[l] = newname
		leaf_dates[newname] = date
		ll.name = newname
	return node
##
#rename leaves
tree = rename_leaves(tree)
	
######Russian leaves
ru_leaves = list()
for leaf in tree.iter_leaf_names():
	if ("Russia" in leaf) or ("FMBA" in leaf):
		ru_leaves.append(leaf)
#print ru_leaves

##########################
### remove all branches that are shorter than the threshold value

def convert_to_nonbinary(t, threshold):
	for node in t.iter_descendants("postorder"):
		#print node.dist
		if node.dist < threshold:
			if not node.is_leaf():
				for child in node.children:
					(node.up).add_child(child, dist = child.dist)
				node.detach()
	return(t)

if nbchecker == 0: #check whether the input tree is nonbinary
	tree = convert_to_nonbinary(tree, cutoff)
	tree.write(format=1, outfile=options.infile+"_nonbinary.nwk")

##########################
##get clusters
clusters = {}

num = 1

#print tree.children

for node in tree.iter_descendants("postorder"):
	if node.is_leaf():
		if node.name in ru_leaves:
			node.add_feature('state',1)
		else:
			node.add_feature('state',0)
		
		#print node.name, node.state
	else:
		node.state = 1
		rus = []
		nonrus = 0
		for c in node.children:
			if c.state == 1:
				rus.append(c)
			else:
				nonrus +=1
		if nonrus > 0:
			node.state = 0
			if len(rus) > 0:
				for child in rus:
					if not child.is_leaf() and child.dist > cutoff:   ### the child is a cluster founder if it is internal (has at least two terminal descendants) and russian
						clusters[child.name] = child.get_leaf_names()
						num +=1

############################################################################################
###################################

##
def addfeatures(node):
	for f in node:
		if "/" in f.name: 
			country = f.name.split("/")[0]
			if (country.isupper()) & ("USA" not in country) & ("DRC" not in country):
				country = country.title()
			if "_" in country:
				country = country.replace("_", " ")
			if "?" in country:
				country = ""
			f.add_feature("country",country)
	return node	
##
def get_earliest_date(node):
	predate = datetime.strptime("2300-01-01", "%Y-%m-%d")
	for l in node:
		if leaf_dates[l.name] < predate:
			predate = leaf_dates[l.name]
	return predate
##############################
##############################
#Write clusters to a file
with open(options.infile+".clusters", 'w') as outf:
	for k in clusters:
			####print clusters
		for elem in clusters[k]:
			if options.verbose:
				print (str(k)+"\t"+elem.replace("\'",""))
			outf.write(str(k)+"\t"+elem.replace("\'","")+"\n")	
##############################
#Visualize one cluster
##############
#manually input cluster and color
clid = ""
if options.cluster == '':
	strclid = input("please enter id:")
	clid = strclid
else:
	clid = options.cluster
##############
common_node = tree.get_common_ancestor(clusters[clid])
#print(common_node)
clusterDate = get_earliest_date(common_node)
########
if options.verbose:
	print(clid+"\t"+str(clusterDate.strftime("%Y-%m-%d")))
########
if options.verbose:
	print(clusterDate)
#common_node.ladderize()
up = common_node.up
#up.show()
up = addfeatures(up)

minimal_dist = common_node.dist

#print(minimal_dist)

#get set of countries on the stem and further
leaves_on_stem_anc = list()
leaves_on_stem_not = list()
distant_leaves = list()
country_min_date = dict()
for leaf in up:
	if leaf not in common_node.get_leaves():
		dist = common_node.get_distance(leaf.name)
		if dist <= (minimal_dist):
			if (leaf_dates[leaf.name] < clusterDate)  or ("Russia" in leaf.name):
				leaves_on_stem_anc.append(leaf.country)
				if (not leaf.country in country_min_date.keys()) or (country_min_date[leaf.country] > leaf_dates[leaf.name]): 
					country_min_date[leaf.country] =leaf_dates[leaf.name]
				#print(str(leaf_dates[leaf.name]))
			else:
				leaves_on_stem_not.append(leaf.country)
		else:
			distant_leaves.append(leaf.country)


##
def add_one(li):
	if len(li) == 0:
		li.append(" ")
	return li
###set up regions
Europe = ["Austria", "Belarus","Moldova","Belgium", "Bosnia and Herzegovina","Bulgaria","Croatia", "Cyprus","Czech Republic", "Denmark", "Finland", "France", "Estonia","Germany", "Scotland","Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Montenegro","Netherlands", "north Macedonia", "Norway", "Poland", "Portugal", "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", "Turkey", "Northern Ireland","Wales","England", "Ukraine"]
NorthAmerica = ["Canada", "Puerto Rico","Costa Rica", "Jamaica", "Mexico", "Panama", "USA"]
Asia = ["Bangladesh","Brunei","Cambodia","China","Georgia","Hong Kong","India","Indonesia","Iran","Israel","Japan","Jordan","Kazakhstan","Kuwait","Lebanon","Malaysia","Nepal","Oman","Pakistan","Philippines","Saudi Arabia","SaudiArabia","Singapore","South Korea","Sri Lanka","Taiwan","Thailand","Timor-Leste","United Arab Emirates","Vietnam", "Abu Dhabi", "Africa", "Aichi", "Amman", "Andhra Pradesh", "Anhui", "Ashdod", "Assam", "Bahrain", "Bangkok", "Bangladesh", "Banten", "Barishal", "Bat Yam", "Beijing", "NanChang", "Hangzhou", "Beirut", "Beit Shemesh", "Beni Brak", "Bihar", "Brunei", "Buraimi", "Capital Governorate", "Central District", "Central Java", "Chattogram", "Chiba", "Chongqing", "Chubu", "Chungcheongnam", "Chushikoku", "Colombo", "Dakhiliyah", "Delhi", "Dhahirah", "Dhaka", "Dubai", "East Java", "East Kalimantan", "Elad", "Elkana", "Fujian", "Gani Tikva", "Gedera", "Gilgit Baltistan", "Guangdong", "Gujarat", "Gyeonggi Province", "Haifa", "Haifa District", "Haryana", "Heilongjiang", "Henan", "Ho Chi Minh City", "Hokkaido", "Hokuriku", "Hokurikushinsyu", "Holon", "Hong Kong", "Hubei", "Hukuk", "India", "Iran", "Irbid", "Ishikawa", "Islamabad", "Israel", "Jakarta", "Jammu", "Japan", "Jeddah", "Jerusalem District", "Jiangsu", "Jiangxi", "Kalutara", "Kanagawa", "Kansai", "Kanto", "Karachi", "Karnataka", "Kathmandu", "Keelung", "Kfar Habad", "Khulna", "Kiryat Gat", "Kiryat Ono", "Kochi", "Kuala Lumpur", "Kuwait", "Kyoto", "Kyusyu", "Ladakh", "Lao Cai", "Liaoning", "Maalot Tarshicha", "Madhya Pradesh", "Madinah", "Maharashtra", "Makkah", "Malaysia", "Manila", "Mongolia", "Muscat", "Mymensingh", "Nagasaki", "Nakhonnayok", "Nara", "National Capital Region", "Netanya", "Netivot", "New Taipei City", "Nonthaburi", "North Batinah", "North District", "North Sharqiyah", "North Sulawesi", "Nur-Sultan", "Odisha", "Osaka", "Pahang", "Pakistan", "Pardesia", "Pathum Thani", "Petah Tikva", "Philippines", "Phuket", "Punjab", "Quangninh", "Raanana", "Rajasthan", "Rajshahi", "Ramat Gan", "Ramla", "Rangpur", "Rawalpindi", "Red River Delta", "Riyadh", "Saitama", "Samut Prakarn", "Saudi Arabia", "Selangor", "Semnan", "Seoul", "Shandong", "Shanghai", "Sharja", "Shoam", "Sichuan", "Sihanoukville", "Singapore", "South Batinah", "South Coast District", "South District", "South Korea", "Special Region of Yogyakarta", "Sri Lanka", "Sylhet", "Tainan", "Taiwan", "Tamil Nadu", "Tbilisi", "Tehran", "Tel Aviv District", "Telangana", "Thailand", "Thanh Hoa", "Timor-Leste", "Tirat Zvi", "Tohoku", "Tokyo", "Trang", "Tveria", "United Arab Emirates", "Uttar Pradesh", "Uttarakhand", "Vietnam", "Vinhphuc", "West Bengal", "West Java", "Yangon", "Yavne", "Yehud", "Yunnan", "Zefat", "Zhejiang"]
SouthAmerica = ["Argentina","Brazil","Chile","Colombia","Ecuador","Peru","Uruguay", "Suriname","Venezuela"]
Africa = ["Algeria", "Asia", "Benin", "Democratic Republic of the Congo", "DRC", "Egypt", "Gambia", "Ghana", "Kenya", "Madagascar", "Mali", "Morocco", "Nigeria", "Senegal", "Sierra Leone", "South Africa", "Tunisia", "Uganda", "Zambia"]
Oceania = ["Auckland", "Australia", "Bay of Plenty", "Canterbury", "Combined Wellington", "Counties Manukau", "Guam", "Lakes", "Midcentral", "Nelson Marlborough", "New South Wales", "New Zealand", "Northern Territory", "Otago", "Queensland", "South Australia", "Southern", "Tasmania", "Victoria", "Waikato", "Wairarapa", "Waitemata", "Wellington", "Western Australia"]
reg_dict = dict()
for c in Europe:
	reg_dict[c] = "Europe"
for c in NorthAmerica:
	reg_dict[c] = "North America"
for c in Asia:
	reg_dict[c] = "Asia"
for c in SouthAmerica:
	reg_dict[c] = "South America"
for c in Africa:
	reg_dict[c] = "Africa"
for c in Oceania:
	reg_dict[c] = "Oceania"
##
def set_region(li,num):
	reli = list()
	if len(set(li)) >num:
		for el in li:
			if el in reg_dict:
				reli.append(reg_dict[el])
			else:
				reli.append(el)
		return reli
	else:
		return li
###
def merge_regions(li, num):
	returnli = list()
	if len(set(li)) >= num:
		returnli = ['Rest of the world']
	else:
		returnli = li
	return returnli
		
###
def set_anc_countries_dates_and_merge(li,num):
	newdict = dict()
	returnli =list()
	if len(set(li)) >= num:
		predate = datetime.strptime("2300-01-01", "%Y-%m-%d")
		for el in li:
			if el in reg_dict:
				if (reg_dict[el] not in newdict.keys()) or (country_min_date[el] < predate):
					newdict[reg_dict[el]] = country_min_date[el]
					predate = country_min_date[el]
			else:
				newdict[el] = country_min_date[el]
		for k in newdict:
			returnli.append(k+"|"+newdict[k].strftime("%b-%d"))
	else:
		for el in set(li):
			if el != ' ':
				returnli.append(el+"|"+str(country_min_date[el].strftime("%b-%d")))
			else:
				returnli.append(" ")
	return returnli
	
	##
	for el in list(set(reli)):
		if el in reg_dict.keys():
			predate = datetime.strptime("2300-01-01", "%Y-%m-%d")
			for countries in reg_dict[el]:
				if countries in country_min_date.keys():
					if country_min_date[countries] < predate:
						predate = country_min_date[countries]
			returnli.append(el+"|"+str(predate.strftime("%b-%d")))
		else:
			#print("IN")
			if el != ' ':
				returnli.append(el+"|"+str(country_min_date[el].strftime("%b-%d")))
			else:
				returnli.append(" ")
				
	return returnli

##
leaves_on_stem_anc = add_one(leaves_on_stem_anc)
leaves_on_stem_anc_corr_date = set_anc_countries_dates_and_merge(leaves_on_stem_anc, 8)
###########
if "Russia" in country_min_date:
	print(clid+"\t"+str(country_min_date["Russia"].strftime("%Y-%m-%d")))
else:
	print(clid+"\tNo Russia on stem")
#########
#leaves_on_stem_not = add_one(leaves_on_stem_not)
distant_leaves = add_one(distant_leaves)
#leaves_on_stem_not = set_region(leaves_on_stem_not,2)
#leaves_on_stem_not = merge_regions(leaves_on_stem_not,2)
distant_leaves = set_region(distant_leaves,2)

######create new tree
dist_leaf_name = "\n".join(list(set(distant_leaves)))
close_leaf_name = "\n".join(list(set(leaves_on_stem_anc_corr_date)))+"\n"
if len(leaves_on_stem_not) >0:
	close_leaf_name += "Non-Russian sequences with later dates ("+str(len(leaves_on_stem_not))+")"
#print (list(set(leaves_on_stem_not)))

##new tree
newt = ete3.Tree(format=1)
RO = newt.add_child(dist = 1.0)
CLL = RO.add_child(name = close_leaf_name,dist=0.0)
DL = RO.add_child(name = dist_leaf_name,dist=0.5)
R = RO.add_child(common_node, name = "R", dist = minimal_dist)

if options.verbose:
	print(newt)

for le in newt:
	if "\n" in le.name:
		newname = le.name.replace("\n","XXXX")
		le.name = newname

newt.write(format=1, outfile=clid+".nwk")


