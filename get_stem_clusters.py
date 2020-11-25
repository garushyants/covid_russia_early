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
parser.add_option('-l', '--leafnames', default= '',help='', type='str')
parser.add_option('-t', '--nonbinary', default= 1, help='', type='int')
parser.add_option('-r', '--root', default= "hCoV-19_Wuhan-Hu-1_2019_EPI_ISL_402125_2019-12-31", help='', type='str')
#parser.add_option('-n','--number',help='',type='int')
parser.add_option('-b','--bootstrap',default = 1, type = 'int')

##get options
options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format=options.bootstrap)
root = options.root
nbchecker = options.nonbinary
#clid=options.number
if options.leafnames == '':
	print("No file with leaves names provided")	
outfile = options.infile+".unique_singletons.csv"
###
cutoff = 1e-05

####
#reroot
tree.set_outgroup(tree&root)
####
#Read correct names
leaf_dates = dict()
leaf_names = dict()
with open(options.leafnames, 'r') as leafnames:
	for l in leafnames:
		l = l[:-1]
		parts = l.split(" -> ")
		pattern = re.compile("hCoV-19/", re.IGNORECASE)
		prename = pattern.sub("", parts[0])
		#prename = parts[0].replace("hCoV-19/","")
		name,gisid,datestr = prename.split("|")
		datecheck = datestr.split("-")
		if len(datecheck)<3:
			if len(datecheck) ==2:
				datestr = datestr+"-28"#set to max possible date
			if len(datecheck) ==1:
				datestr = datestr+"-12-31"#set to the latest possible date if data is missing
		date = datetime.strptime(datestr, "%Y-%m-%d")
		newdateformat = date.strftime("%b-%d")
		newname=""
		if "Russia" in l:
			newname = name+"|"+str(newdateformat)
		else:
			newname = name+"|"+gisid+"|"+str(newdateformat)
		leaf_names[parts[1]] = newname
		leaf_dates[newname] = date
##
def rename_leaves(node):
	for l in node:
		l.name = leaf_names[l.name]
	return node
##
#rename leaves
tree = rename_leaves(tree)
	
######Russian leaves
ru_leaves = list()
for leaf in tree.iter_leaf_names():
	if "Russia" in leaf:
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
			if "Russia" in l.name and (l.name !=intname and l.up == ancestor):
				boolcheck = False
				break
			elif "Russia" not in l.name:
				intcheck = True
		if boolcheck and intcheck:
			int_cases.append(intname)

int_cases.sort()
###Print to file	
with open(outfile, 'w') as outf:	
	counter = 1
	for el in int_cases:
		outf.write("UniqueSingletons"+str(counter)+"\t"+el+"\n")
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
def layout(node):
	if node.is_leaf():
		if "Russia/" in node.name:
			Rudesc = faces.AttrFace("name", fsize=10)
			Rudesc.margin_left = 100
			Rudesc.background.color = '#bdbdbd'
			faces.add_face_to_node(Rudesc, node, 0)
		else:
			if "_" in node.name:
				node.name = node.name.replace("_"," ")
			countrydesc = faces.AttrFace("name", fsize=10)
			countrydesc.margin_left = 5
			countrydesc.margin_bottom = 10
			countrydesc.margin_top = 5
			countrydesc.margin_right = 5
			faces.add_face_to_node(countrydesc, node, 0, aligned = False)
##
def addfeatures(node):
	for f in node:
		if "/" in f.name: 
			country = f.name.split("/")[0]
			f.add_feature("country",country)
	return node	

#########################################
#Visualize one cluster
###Select cluster from keyboard
strclid = input("please enter id:")
clid = int(strclid) - 1

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
Europe = ["Austria", "Belgium", "Croatia", "Czech_Republic", "Denmark", "Finland", "France", "Estonia","Germany", "Scotland","Greece", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Netherlands", "Norway", "Poland", "Portugal", "Romania", "Serbia", "Slovakia", "Slovenia", "Spain", "Sweden", "Switzerland", "Turkey", "Northern_Ireland","Wales","England"]
NorthAmerica = ["Canada", "Puerto_Rico","Costa_Rica", "Jamaica", "Mexico", "Panama", "USA"]
Asia = ["Bangladesh","Brunei","Cambodia","China","Georgia","Hong_Kong","India","Indonesia","Iran","Israel","Japan","Jordan","Kazakhstan","Kuwait","Lebanon","Malaysia","Nepal","Oman","Pakistan","Philippines","Saudi_Arabia","Singapore","South_Korea","Sri_Lanka","Taiwan","Thailand","Timor-Leste","United_Arab_Emirates","Vietnam"]
SouthAmerica = ["Argentina","Brazil","Chile","Colombia","Ecuador","Oceania","Peru","Uruguay"]
Africa = ["DRC","Egypt","Gambia","Ghana","Kenya","Morocco","Nigeria","Senegal","South_Africa","Tunisia","Uganda"]
reg_dict = dict()
for c in Europe:
	reg_dict[c] = "Europe"
for c in NorthAmerica:
	reg_dict[c] = "North_America"
for c in Asia:
	reg_dict[c] = "Asia"
for c in SouthAmerica:
	reg_dict[c] = "South_America"
for c in Africa:
	reg_dict[c] = "Africa"
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
	close_leaf_name += "Sequences with later dates ("+str(len(leaves_on_stem_not))+")"

newt = ete3.Tree(format=0)
R = newt.add_child(name = "R", dist = up.dist)#, dist = int_leaf.get_distance(int_leaf.name))
CLL = R.add_child(name = close_leaf_name, dist=0.0)
#print(newt)
DL = R.add_child(name = dist_leaf_name, dist=0.5)
INT = R.add_child(name = int_leaf.name, dist = minimal_dist)
#newt.set_outgroup(R)
print(newt)
#newt.show()

##visualize
ts = TreeStyle()
ts.branch_vertical_margin = 1
ts.root_opening_factor = 1
ts.scale = 50
ts.show_leaf_name = False
ts.layout_fn = layout
#nstyle = NodeStyle()
rnstyle = NodeStyle()
lnstyle = NodeStyle()
dlstyle = NodeStyle()
#color = random_color()
for rn in newt.traverse():
	rnstyle["vt_line_width"] = 1
	rnstyle["hz_line_width"] = 1
	rnstyle["size"] = 0
	rn.set_style(rnstyle)
	if rn.is_leaf() and "Russia" in rn.name:
		#lnstyle["bgcolor"] = '#bdbdbd'
		lnstyle["hz_line_width"] = 1
		lnstyle["size"] = 5
		lnstyle["fgcolor"] = "black"
		rn.set_style(lnstyle)
###set distant leaves style
dlstyle["hz_line_width"] = (DL.name.count("\n")+1)*15
dlstyle["size"] = 0
DL.set_style(dlstyle)
#newt.show()
newt.render("Figures/Singletons_unique_"+str(clid+1)+".pdf", w = 4000, units= 'px', dpi = 350, tree_style = ts)

