./RigidMutualInfoRegistration Data/HighResCT.hdr Data/LowResCT.hdr Temp/LowResCTReg.hdr Temp/rigidTrans.txt 
./RigidTransformation Data/LowResPET.hdr try.hdr rigidTrans.txt Data/HighResCT.hdr 

./ResampleRatio Data/HighResCT.hdr Temp/Fixed.hdr 1
./ResampleRatio Temp/LowResCTReg.hdr Temp/Moving.hdr 1

./MSSOFuzzyConnSeg3DInt Data/seed100.txt Temp/movingFC.hdr Temp/Moving.hdr 100 -900 3
./MSSOFuzzyConnSeg3DInt Data/seed100.txt Temp/fixedFC.hdr Temp/Fixed.hdr 0 -900 3

./Threshold Temp/movingFC.hdr Temp/movingBin.hdr 1000 100000
./Threshold Temp/movingFC.hdr Temp/movingBinT.hdr 9000 100000

./Threshold Temp/fixedFC.hdr Temp/fixedBin.hdr 1000 100000
./Threshold Temp/fixedFC.hdr Temp/fixedBinT.hdr 9700 100000

./BinaryDilate Temp/movingBinT.hdr Temp/movingBinT.hdr 3
./MaskInverse Temp/movingBinT.hdr Temp/movingBinT.hdr
./MaskIntersection Temp/movingBin.hdr Temp/movingBinT.hdr Temp/movingBin.hdr
./ConnectComponent Temp/movingBin.hdr Temp/movingBinCC.hdr 1000 10000000000
./Threshold Temp/movingBinCC.hdr Temp/movingBin.hdr 1 1
./BinaryDilate Temp/movingBin.hdr Temp/movingBin.hdr 30
./BinaryErode Temp/movingBin.hdr Temp/movingBin.hdr 30

./BinaryDilate Temp/fixedBinT.hdr Temp/fixedBinT.hdr 1
./MaskInverse Temp/fixedBinT.hdr Temp/fixedBinT.hdr
./MaskIntersection Temp/fixedBin.hdr Temp/fixedBinT.hdr Temp/fixedBin.hdr
./ConnectComponent Temp/fixedBin.hdr Temp/fixedBinCC.hdr 1000 10000000000
./Threshold Temp/fixedBinCC.hdr Temp/fixedBin.hdr 1 1
./BinaryDilate Temp/fixedBin.hdr Temp/fixedBin.hdr 20
./BinaryErode Temp/fixedBin.hdr Temp/fixedBin.hdr 20

./MaskUnion Temp/movingBin.hdr Temp/fixedBin.hdr Temp/whole.hdr
./BinaryDilate Temp/whole.hdr Temp/whole.hdr 5

./ExtractMaskROI Temp/Moving.hdr Temp/whole.hdr Temp/movingImg.hdr
./ExtractMaskROI Temp/movingBin.hdr Temp/whole.hdr Temp/movingMsk.hdr

./ExtractMaskROI Temp/Fixed.hdr Temp/whole.hdr Temp/fixedImg.hdr
./ExtractMaskROI Temp/fixedBin.hdr Temp/whole.hdr Temp/fixedMsk.hdr

./MultiScaleHessianToObject Temp/movingImg.hdr Temp/movingVes.hdr Temp/movingSca.hdr 1 0.5 2.5 5 1
./MultiScaleHessianToObject Temp/fixedImg.hdr Temp/fixedVes.hdr Temp/fixedSca.hdr 1 0.5 2.5 5 1

./LungDeformReg parameters.txt

./ReferenceRewrite Result/Deformed12_o_1_1_1_4_4_4.mha Temp/fixedImg.hdr Temp/regImg.hdr

./imageDeform Temp/movingImg.hdr Result/Disp12_x_1_1_1_4_4_4.mha Result/Disp12_y_1_1_1_4_4_4.mha Result/Disp12_z_1_1_1_4_4_4.mha Temp/movingWrp.mha

./ReferenceRewrite Temp/movingWrp.mha Temp/fixedImg.hdr Temp/movingWrp.hdr
