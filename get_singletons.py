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
parser.add_option('-c','--cluster', default = -2, help='',type='int')
parser.add_option('-b','--bootstrap',default = 1, type = 'int')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.bootstrap)
#root = options.root
nbchecker = options.nonbinary
#clid=options.number
outfile = options.infile+".singletons.csv"
###
cutoff = 1e-05

####
#reroot
#tree.set_outgroup(tree&root)
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
	if ("Russia" in leaf) or ("FMBA" in leaf): #Some FMBA sequences does not have country in metadata
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
	#tree.write(format=1, outfile=options.infile+"_nonbinary.nwk")


##########################
## get interesting cases
int_cases = list()
for leaves in tree:
	intname = leaves.name
	if intname in ru_leaves:
		boolcheck = True
		intcheck = False
		ancestor = leaves.up
		for l in ancestor:
			if ("Russia" in l.name or "FMBA" in l.name): 
				if (l.name !=intname and l.up == ancestor):
					boolcheck = False
					break
			else:
				intcheck = True
		if boolcheck and intcheck:
			int_cases.append(intname)

int_cases.sort()
###Print to file	
with open(outfile, 'w') as outf:	
	counter = 1
	for el in int_cases:
		outf.write("Singletons"+str(counter)+"\t"+el+"\n")
		counter += 1		
####
print("\n".join(int_cases))
						 
############################################################################################
###################################
#####Visualization
###################################
def random_color():
    color = "#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
    return color

####
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

#########################################
#Visualize one cluster
###Select cluster from keyboard
clid = -1
if options.cluster == -2:
	strclid = input("please enter id:")
	clid = int(strclid) -1
else:
	clid = options.cluster -1
##########
int_leaf = tree&int_cases[clid]
clusterDate = leaf_dates[int_cases[clid]]
print(clusterDate)
print(int_leaf)
up = int_leaf.up
up = addfeatures(up)

minimal_dist = int_leaf.dist

#get set of countries on the stem and further
leaves_on_stem_anc = list()
leaves_on_stem_not = list()
distant_leaves = list()
country_min_date = dict()
for leaf in up:
	if leaf != int_leaf:
		dist = int_leaf.get_distance(leaf.name)
		if dist <= (minimal_dist + cutoff):
			if leaf_dates[leaf.name] < clusterDate:
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
			print("IN")
			if el != ' ':
				returnli.append(el+"|"+str(country_min_date[el].strftime("%b-%d")))
			else:
				returnli.append(" ")
				
	return returnli
##
leaves_on_stem_anc = add_one(leaves_on_stem_anc)
leaves_on_stem_anc_corr_date = set_anc_countries_dates_and_merge(leaves_on_stem_anc, 8)
#leaves_on_stem_not = add_one(leaves_on_stem_not)
distant_leaves = add_one(distant_leaves)
#leaves_on_stem_not = set_region(leaves_on_stem_not,2)
#leaves_on_stem_not = merge_regions(leaves_on_stem_not,2)
distant_leaves = set_region(distant_leaves,2)
#print (set(leaves_on_stem))
#print (set(distant_leaves))

######create new tree
dist_leaf_name = "\n".join(list(set(distant_leaves)))
close_leaf_name = "\n".join(list(set(leaves_on_stem_anc_corr_date)))+"\n"
if len(leaves_on_stem_not) >0:
	close_leaf_name += "Non-Russian sequences with later dates ("+str(len(leaves_on_stem_not))+")"

##
newt = ete3.Tree(format=1)
RO = newt.add_child(dist = up.dist)
CLL = RO.add_child(name = close_leaf_name,dist=0.0)
DL = RO.add_child(name = dist_leaf_name,dist=0.5)
R = RO.add_child(name = int_leaf.name, dist = minimal_dist)

print(newt)

for le in newt:
	if "\n" in le.name:
		newname = le.name.replace("\n","XXXX")
		le.name = newname

newt.write(format=1, outfile="Singleton"+str(clid+1)+".nwk")
