#!/usr/bin/python3

import sys
import ete3
import optparse
from ete3 import Tree, NodeStyle, AttrFace, faces, TreeStyle, SeqMotifFace

parser=optparse.OptionParser()
parser.add_option('-i', '--infile', help='Input tree file', type='str')
parser.add_option('-f','--outfolder',default = "Figures/", type = 'str')

options, args=parser.parse_args()
tree = ete3.Tree(options.infile, format =1) #input tree is from treetime and nonbinary

###
for l in tree:
	if "XXXX" in l.name:
		newname = l.name.replace("XXXX", "\n")
		l.name = newname
###

print(tree)

DL = tree.search_nodes(dist = 0.5)[0]


def layout(node):
	if node.is_leaf():
		if ("Russia/" in node.name) or ("FMBA" in node.name):
			Rudesc = faces.AttrFace("name", fsize=10)
			Rudesc.margin_left = 5
			Rudesc.margin_right = 5
			Rudesc.background.color = "#d9d9d9"
			faces.add_face_to_node(Rudesc, node, 0, aligned=True)
		else:
			if "_" in node.name:
				node.name = node.name.replace("_"," ")
			countrydesc = faces.AttrFace("name", fsize=10)
			countrydesc.margin_left = 5
			countrydesc.margin_bottom = 10
			countrydesc.margin_top = 5
			countrydesc.margin_right = 5
			faces.add_face_to_node(countrydesc, node, 0, aligned = False)

#########################
##visualize
##visualize
ts = TreeStyle()
ts.branch_vertical_margin = 1
ts.root_opening_factor = 1
ts.scale = 50
ts.draw_guiding_lines =True
ts.guiding_lines_type = 2

ts.show_leaf_name = False
ts.layout_fn = layout
#nstyle = NodeStyle()
rnstyle = NodeStyle()
lnstyle = NodeStyle()
dlstyle = NodeStyle()
#color = random_color()
for rn in tree.traverse():
	rnstyle["vt_line_width"] = 1
	rnstyle["hz_line_width"] = 1
	rnstyle["size"] = 0
	rn.set_style(rnstyle)
	if rn.is_leaf() and "Russia/" in rn.name:
		#lnstyle["bgcolor"] = '#bdbdbd'
		lnstyle["hz_line_width"] = 1
		lnstyle["size"] = 5
		lnstyle["fgcolor"] = "black"
		rn.set_style(lnstyle)
###set distant leaves style
dlstyle["hz_line_width"] = (DL.name.count("\n")+1)*15
dlstyle["size"] = 0
DL.set_style(dlstyle)

temp = options.infile
outfile = temp.replace(".nwk",".pdf")
tree.render(options.outfolder+"/"+outfile, w = 4000, units= 'px', dpi = 350, tree_style = ts)

