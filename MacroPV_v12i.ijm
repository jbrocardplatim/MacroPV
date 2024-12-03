//Macro used to measure co-labelling of WFA, tomato and PV in cells from brain slices
//Needs the Cellpose wrapper installed : https://github.com/BIOP/ijl-utilities-wrappers
//@Jacques Brocard for Juan Ren & Sabine Richard (IGFL), 2024


//---INITIALIZATION
run("Set Measurements...", "area mean redirect=None decimal=3");
//Microscope where the images have been taken => will change expected intensities
Dialog.create("Images obtained...");
Dialog.addChoice("Microscope",newArray("LSM780", "LSM800"));
Dialog.show();
micro = Dialog.getChoice();

cell_radius=15; //Mean cell radius in pixels
nb_z=5; //Nb of planes
n_c1="WFA"; //Labelling to be detected in 1st channel
n_c2="tomato"; //Labelling to be detected in 2nd channel
n_c3="PV"; //Labelling to be detected in 3rd channel

//---PICK DIRECTORY and analyzes each .czi stack
dir = getDirectory("Choose the Directory");
list = getFileList(dir);
if (list.length>0) print("Microscope\t Filename\t #cells \ttotal "+ n_c1 +" \ttotal "+ n_c2 +"\ttotal "+ n_c3 +"\ttotal "+ n_c1+"+"+n_c2+" \ttotal "+ n_c1+"+"+n_c3+" \ttotal "+n_c2+"+"+n_c3+" \tAll labelings \tArea (µm²) \traw WFA Signal");
for (i=0; i<list.length; i++){
    ext=substring(list[i],lengthOf(list[i])-4,lengthOf(list[i]));
    if (ext==".czi"){
    	run("Bio-Formats Windowless Importer", "open=["+dir+list[i]+"]");
    	getDimensions(w, h, channels, slices, frames);
    	nb_c=channels;
    	roiManager("Reset");
    	main(list[i]);
    	close();
    }
    selectWindow("Log");
	saveAs("Text",dir+"Results_v12i.txt");
} 


function main(t){

	//--- INITIALIZE STACK
	mean=newArray(nb_z);
	t=getTitle();
	dir=getDirectory("image");
		
	stack=substring(t,0,lengthOf(t)-4)+"_stack.tif";
	temp=substring(t,0,lengthOf(t)-4)+"_temp.tif";
	
	//Automatic projection of the three most intense z planes
	for (s=0;s<floor(nSlices/nb_c);s++){
		setSlice(nb_c+s*nb_c);
		getStatistics(area, mean[s]);
	}
	maxi=0;
	for (s=1;s<nb_z;s++){
		if (mean[s]>mean[maxi]) maxi=s;
	}
	deb=maxi-1;
	fin=maxi+1;
	if (deb<0) {
		deb++; fin++;
	}
	if (fin>nb_z) {
		deb--; fin--;
	}
	run("Z Project...", "start="+deb+" stop="+fin+" projection=[Max Intensity]");
	
	//Images obtained with LSM800 need to be adjusted
	if (micro=="LSM800"){
		run("Scale...", "x=0.75 y=0.75 z=1.0 interpolation=Bilinear average create");
		run("Arrange Channels...", "new=132");	
	}
	//Save reference original stack
	rename(stack);
	saveAs("Tiff",dir+stack);
	
	setSlice(1);
	//In case there is a control DAPI image, erase it
	if (nb_c==4) run("Delete Slice", "delete=channel");
	
	//Get original area and mean WFA signal for later normalization
	setSlice(1);
	run("Select All");
	getStatistics(Area_WFA,Mean_WFA);
	run("Properties...", "channels=3 slices=1 frames=1 pixel_width=1 pixel_height=1 voxel_depth=1 unit=pixel");
	//Smooth and subtract background before image analysis
	run("Smooth","stack");
	if (micro=="LSM800") {
		run("Subtract Background...", "rolling="+2*cell_radius+" stack");
	}else {
		run("Subtract Background...", "rolling="+cell_radius+" stack");
	}
	//Save analysis-ready stack as a "_temp.tif" file
	saveAs("Tiff",dir+temp);

	//Automatic detection of cell contours using Cellpose with channel 2 = "tom" labelling
	//and creates a mask with the extract_cellpose function
	selectWindow(temp);
	run("Cellpose Advanced", "diameter="+2*cell_radius+" cellproba_threshold=0.5 flow_threshold=0.4 anisotropy=1.0 diam_threshold=10.0 model=cyto2 nuclei_channel=0 cyto_channel=2 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
	rename("tom");
	extract_cellpose("tom");
	
	//Automatic detection of cell contours using Cellpose with channel 1 = "wfa" labelling
	//and creates a mask with the extract_cellpose function
	selectWindow(temp);
	run("Cellpose Advanced", "diameter="+2*cell_radius+" cellproba_threshold=0.5 flow_threshold=0.4 anisotropy=1.0 diam_threshold=10.0 model=cyto2 nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
	rename("wfa1");
	extract_cellpose("wfa1");

	//Automatic detection of cell contours using Cellpose with channel2-assisted="tom" channel 1 = "wfa" labelling
	//and creates a mask with the extract_cellpose function
	selectWindow(temp);
	run("Cellpose Advanced", "diameter="+2*cell_radius+" cellproba_threshold=0.5 flow_threshold=0.3 anisotropy=1.0 diam_threshold=10.0 model=cyto2 nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
	rename("wfa2");
	extract_cellpose("wfa2");

	//Automatic detection of cell contours using Cellpose with channel 3 = "pv" labelling
	//and creates a mask with the extract_cellpose function
	selectWindow(temp);
	run("Cellpose Advanced", "diameter="+2*cell_radius+" cellproba_threshold=0.5 flow_threshold=0.2 anisotropy=1.0 diam_threshold=10.0 model=cyto2 nuclei_channel=0 cyto_channel=3 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
	rename("pv");
	extract_cellpose("pv");
		
	//Add up all the masks, separate individual cells and analyze particles
	imageCalculator("Add create", "tom.tif","pv.tif");
	rename("tom_pv");
	imageCalculator("Add create", "wfa1.tif","wfa2.tif");
	rename("wfa");
	imageCalculator("Add create", "tom_pv","wfa");
	run("Watershed");
	cell_size=PI*cell_radius*cell_radius/4;
	run("Analyze Particles...", "size="+cell_size+"-Infinity clear add exclude");
	
	close(); close("tom_pv"); close("wfa");
	close("tom.tif"); close("wfa1.tif"); close("wfa2.tif"); close("pv.tif");

	//Saturate images to similar levels...
	selectWindow(temp);
	//tem_c variables keep the intensity used to saturate said pixels
	temp_c1=saturate_slice(1,0.001);
	temp_c2=saturate_slice(2,0.001);
	temp_c3=saturate_slice(3,0.01);

	//... and classify cells based on individual labelling intensities
	classify_cellpose(temp_c1,temp_c2,temp_c3);
	close(); if (micro=="LSM800") close();
}


function extract_cellpose(te){
	//Create mask from automatically detected cells
	selectWindow(te);
	setAutoThreshold("Default dark");
	setThreshold(1.0000, 1000000000000000000000000000000.0000);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Watershed");
	er=1;
	//Erode differently the cell contours, depending on the channel used for detection
	if (te=="pv") er=2;
	if (te=="wfa2") er=4;
	if (te!="tom") run("Options...", "iterations="+er+" count=1 black pad do=Erode");
	run("Options...", "iterations=1 count=1 black pad do=Open");
	saveAs("Tiff", dir+te+".tif");
}


function saturate_slice(s,percent) { 
	//Modify image display to indicated percentage of saturated pixels, before applying LUT
	setSlice(s);
	run("Select All");
	getHistogram(values, counts, 256);
	beg=0; area=0;
	while (beg<=255) {
		area=area+counts[beg];
		beg++;
	}
	beg=1;
	while(counts[beg+1]>counts[beg]) beg++;
	end=255;
	while(sum_counts(end)/area<percent) end--;

	setMinAndMax(values[beg],values[end]);
	run("Apply LUT", "slice");
	return values[end];
}


function sum_counts(e){
	//Count pixels whose intensity is higher than "e"
	sum=0;
	for (f=255;f>=e;f--){
		sum=sum+counts[f];
	}
	return sum;
}


function classify_cellpose(t1,t2,t3){
	
	//At this step, the ROI Manager is full of detected cell contours
	//c1 to c7 = nb of cells of each combination of labelings
	c1=0; c2=0; c3=0; c4=0; c5=0; c6=0; c7=0;
	//"s_c" are variability thresholds (see below)
	s_c1=0.5;
	s_c2=1;
	s_c3=0.5;
	//Thresholds are slightly higher with LSM800 images since background substraction was not as sharp (lines 86-92)
	th_c1=10000; 
	if (micro=="LSM800") th_c1=12000;
	th_c2=14000;  
	if (micro=="LSM800") th_c2=18000;
	th_c3=22000; 
	
	selectWindow(temp);
	nROIs=roiManager("Count");
	mean=newArray(4);
	sd=newArray(4);
	//For each ROI = cell contour
	for (r=0;r<nROIs;r++){
		roiManager("Select",0);
		for (s=1;s<=3;s++){
			setSlice(s);
			//Mean and SD of each channel is measured...
			getStatistics(a,mean[s],b,c,sd[s]);
		}
		col="black"; c1pos=false; c2pos=false; c3pos=false;
		
		//... and compared to mean and variability thresholds to establish identity
		if (mean[1]>th_c1) c1pos=true; 
		//c1 = "WFA", must be very variable
		if (mean[1]>th_c1/2 && sd[1]/mean[1] > s_c1) c1pos=true;
		//If channel 1 saturation threshold was too low, discard current ROI
		if (t1<6000) c1pos=false; 
		
		//c2 = "tom" 
		if (mean[2]>th_c2) c2pos=true; 
		if (mean[2]>th_c2/2 && sd[2]/mean[2] < s_c2) c2pos=true; 
		//If channel 2 saturation threshold was too low, discard current ROI
		if (t2<6000) c2pos=false; 
		
		//c2 = "pv" 
		if (mean[3]>th_c3) c3pos=true; 
		if (mean[3]>th_c3/2 && sd[3]/mean[3] < s_c3) c3pos=true; 
		//If channel 3 saturation threshold was too low, discard current ROI
		if (t3<6000) c3pos=false; 
	
		//Color-code cell contours based on labellings combinations
		if (c1pos & !c2pos & !c3pos) {
			c1++; col="green";
		}
		if (!c1pos & c2pos & !c3pos) {
			c2++; col="red";
		}
		if (!c1pos & !c2pos & c3pos) {
			c3++; col="blue";
		}
		if (c1pos & c2pos & !c3pos) {
			c4++; col="orange";
		}
		if (c1pos & !c2pos & c3pos) {
			c5++; col="cyan";
		}
		if (!c1pos & c2pos & c3pos) {
			c6++; col="magenta";
		}
		if (c1pos & c2pos & c3pos) {
			c7++; col="white";
		}
		roiManager("Rename", (r+1)+"-"+col);
		Roi.setStrokeColor(col);
		roiManager("Add");
		roiManager("Delete");
	}
	
	nb_cells=c1+c2+c3+c4+c5+c6+c7;
	c1=c1+c4+c5+c7;
	c2=c2+c4+c6+c7;
	c3=c3+c5+c6+c7;
	c4=c4+c7;
	c5=c5+c7;
	c6=c6+c7;
	
	//Save color-coded ROI Manager and print cell countings of each identity
	roiManager("Save",dir+substring(stack,0,lengthOf(stack)-4)+"_v12i.zip");
	print(micro+" \t"+t+" \t"+nb_cells+" \t"+c1+" \t"+c2+" \t"+c3+" \t"+c4+" \t"+c5+" \t"+c6+" \t"+c7+" \t"+Area_WFA+" \t"+Mean_WFA);

}