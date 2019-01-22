#From binary skeleton, generate tree with different point types (normal-branching-terminal)
./IniSkeleton skeletonOri.hdr skeLabel.hdr
#From 3-type skeleton, generate tree with parent-children relationship and branch label
./GenTree skeLabel.hdr  1
mv tree.txt treeOri.txt
#From original tree, apply dynamic programming to remove small branches, resulting in new binary tree 
./ReconstTree treeOri.txt skeletonOri.hdr skeletonPrune.hdr
#Reconstruct tree to get corrected parent-children relationship and branch label
./IniSkeleton skeletonPrune.hdr skeLabelPrune.hdr  
./GenTree  Temp/skeLabelPrune.hdr 1
mv tree.txt  Temp/treeNew.txt
#Additional function to visualize tree data
#./ImageTree treeOri.txt skeLabel.hdr treeLabelOri.hdr treeParLabelOri.hdr
#./ImageTree treeNew.txt skeLabelPrune.hdr treeLabelNew.hdr treeParLabelNew.hdr