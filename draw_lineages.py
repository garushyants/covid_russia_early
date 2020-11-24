#!/usr/bin/python3

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='Input tree file', type='str')
parser.add_option('-f','--outfolder',default = "Figures_lineages", type = 'str')

options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format =1) #input tree is from treetime and nonbinary

###
for l in tree:
	if "XXXX" in l.name:
		newname = l.name.replace("XXXX", "\n")
		l.name = newname
###

print(tree)

common_node = tree.search_nodes(name= "R")[0]


up = common_node.up

DL = up.children[1]
#print(DL)

def layout(node):
	if node.is_leaf():
		if ("Russia/" in node.name) and ("FMBA" not in node.name):
			Rudesc = faces.AttrFace("name", fsize=10)
			Rudesc.margin_left = 5
			Rudesc.margin_right = 5
			Rudesc.background.color = "#d9d9d9"
			faces.add_face_to_node(Rudesc, node, 0, aligned=True)
		elif ("FMBA" in node.name):
			fmbadesc = faces.AttrFace("name", fsize=10)
			fmbadesc.margin_left = 5
			fmbadesc.margin_right = 5
			fmbadesc.background.color = "#77e8ed"
			faces.add_face_to_node(fmbadesc, node, 0, aligned=True)
		else:
			if "_" in node.name:
				node.name = node.name.replace("_"," ")
			countrydesc = faces.AttrFace("name", fsize=10)
			countrydesc.margin_left = 5
			countrydesc.margin_bottom = 10
			countrydesc.margin_top = 5
			countrydesc.margin_right = 5
			faces.add_face_to_node(countrydesc, node, 0, aligned = False)

color = "#3690c0"#input("please enter color:")
#########################
##visualize
ts = TreeStyle()
ts.branch_vertical_margin = 0.5
ts.root_opening_factor = 1
ts.scale =  50
ts.draw_guiding_lines =True
ts.guiding_lines_type = 2


ts.show_leaf_name = False
ts.layout_fn = layout
nstyle = NodeStyle()
rnstyle = NodeStyle()
lnstyle = NodeStyle()
dlstyle = NodeStyle()
#color = random_color()
for rn in tree.traverse():
	rnstyle["vt_line_width"] = 1
	rnstyle["hz_line_width"] = 1
	rnstyle["size"] = 0
	rn.set_style(rnstyle)
for n in common_node.traverse():
	nstyle["vt_line_color"] = color
	nstyle["hz_line_color"] = color
	nstyle["vt_line_width"] = 4
	nstyle["hz_line_width"] = 4
	nstyle["size"] = 0
	n.set_style(nstyle)
	if n.is_leaf():
		lnstyle["vt_line_color"] = color
		lnstyle["hz_line_color"] = color
		lnstyle["vt_line_width"] = 4
		lnstyle["hz_line_width"] = 4
		#lnstyle["bgcolor"] = "#d9d9d9"
		lnstyle["size"] = 5
		lnstyle["fgcolor"] = "black"
		n.set_style(lnstyle)
###set distant leaves style
dlstyle["hz_line_width"] = (DL.name.count("\n")+1)*15
dlstyle["size"] = 0
DL.set_style(dlstyle)

temp = options.infile
outfile = temp.replace(".nwk",".pdf")
tree.render(options.outfolder+"/"+outfile, w = 4000, units= 'px', dpi = 350, tree_style = ts)

