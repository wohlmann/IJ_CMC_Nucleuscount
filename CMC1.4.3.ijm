macro "n nuclei in elongated cardiomyocytes" {
setBatchMode(true);
dir1 = getDirectory("Please choose source directory ");
list = getFileList(dir1);
dir2 = getDirectory("Please choose destination directory ");
print("\\Clear");
print("Reset: log, Results, ROI Manager");
run("Clear Results");
updateResults;
roiManager("reset");
while (nImages>0) {					
	selectImage(nImages);
	close();
}
print("_");
getDateAndTime(year, month, week, day, hour, min, sec, msec);
print("Starting analysis at: "+day+"/"+month+"/"+year+" :: "+hour+":"+min+":"+sec+"");
print("_");
nImg=(list.length);
IMG=0;
N=0;
for (i=0; i<list.length; i++) {						
	path = dir1+list[i];			
	print("start processing of "+path+"");
	print("_");
	run("Bio-Formats Importer", "open=[path] autoscale color_mode=Default view=Hyperstack stack_order=XYCZT open_all_series");
	title1= getTitle;
	title2 = File.nameWithoutExtension;
	actimage = getTitle();
	N=N+1;
	IMG=IMG+1;
	print("ANALYSING IMAGE "+IMG+" of "+nImg+"");
	print("analysing image "+title1+":");
	print("_");
	run("Duplicate...", "duplicate");
	selectWindow(""+title1+"");
	run("Split Channels");
	selectWindow("C2-"+title1+"");
	run("Duplicate...", " ");
	print("selection of cells:");
	print("deleting overlapping cells");
	selectWindow("C2-"+title1+"");
	setAutoThreshold("IsoData dark");
	run("Convert to Mask");
	run("Analyze Particles...", "size=800-3000 circularity=0.00-0.5 show=Outlines exclude include add");
	roiManager("Select", newArray());
	run("Select All");
	roiManager("Combine");
	roiManager("Add");
	roiManager("deselect");
	ROIc = roiManager("count");
	while (ROIc>1) {
		roiManager("select", 0);
		roiManager("delete");
		ROIc = roiManager("count");
	}
	selectWindow("C2-"+title2+"-1.czi");
	if(ROIc!=0){
		roiManager("select", 0);
		run("Enlarge...", "enlarge=20");
		run("Make Inverse");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
		roiManager("Deselect");
		roiManager("Delete");
		setBackgroundColor(255, 255, 255);
	}
	run("Select None");
	run("Duplicate...", " ");
	print("deleting rounded cells");
	setAutoThreshold("MaxEntropy dark");
	setThreshold(2778, 65535);
	run("Convert to Mask");
	run("Open");
	run("Fill Holes");
	run("Analyze Particles...", "size=300-5000 circularity=0.60-1.00 show=Outlines include add");
	roiManager("Select", newArray());
	run("Select All");
	roiManager("Combine");
	roiManager("Add");
	roiManager("deselect");
	ROIc = roiManager("count");
	while (ROIc>1) {
		roiManager("select", 0);
		roiManager("delete");
		ROIc = roiManager("count");
	}
	selectWindow("C2-"+title2+"-1.czi");
	if(ROIc!=0){
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=10");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
		roiManager("Deselect");
		roiManager("Delete");
		setBackgroundColor(255, 255, 255);
	}
	run("Select None");
	run("Duplicate...", " ");
	print("deleting remaining artifacts");
	setAutoThreshold("Default dark");
	setThreshold(1702, 65535);
	run("Convert to Mask");
	run("Analyze Particles...", "size=10-500 circularity=0.00-1.00 show=Outlines include add");
	roiManager("Select", newArray());
	run("Select All");
	roiManager("Combine");
	roiManager("Add");
	roiManager("deselect");
	ROIc = roiManager("count");
	while (ROIc>1) {
		roiManager("select", 0);
		roiManager("delete");
		ROIc = roiManager("count");
	}
	selectWindow("C2-"+title2+"-1.czi");
	if(ROIc!=0){
		roiManager("Select", 0);
		run("Enlarge...", "enlarge=3");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
		roiManager("Deselect");
		roiManager("Delete");
		setBackgroundColor(255, 255, 255);
	}
	run("Select None");
	run("Duplicate...", " ");
	print("detecting individual cells");
	setAutoThreshold("Default dark");
	setThreshold(1542, 65535);
	run("Convert to Mask");
	run("Analyze Particles...", "size=1000-4000 circularity=0.10-1.00 show=Outlines exclude include add");
	print("fine selection of individual cells");
	ROIc = roiManager("count");
	if(ROIc==0){
		selectWindow("C2-"+title2+"-1.czi");
		run("Select None");
		run("Duplicate...", " ");
		setAutoThreshold("IJ_IsoData dark");
		setThreshold(375, 65535);	
		run("Convert to Mask");
		run("Analyze Particles...", "size=1000-4000 circularity=0.10-1.00 show=Outlines exclude include add");
	}
	print("fine selection of individual cells");
	ROIc = roiManager("count");
	nROI = parseInt(roiManager("count"));
	IJ.renameResults("ResOUT");
	for (x=0; x<nROI; x++) {
		roiManager("Select", x);
		run("Fit Ellipse");
		run("Set Measurements...", "area feret's redirect=None decimal=4");
		run("Measure");
		List.setMeasurements;
		minFer = List.getValue("MinFeret");
		Yfere = List.getValue("FeretY");
		Fer = List.getValue("Feret");
		if ((minFer>=28)|(Yfere >=700)|(Fer <=93)){
			roiManager("Select", x);
			roiManager("Delete");
			nROI = nROI-1;
			x=x-1;
	        }
		}
	print("_");
	print("generating QC image:");
	nROI = parseInt(roiManager("count"));
	for (y=0; y<nROI; y++) {
		roiManager("Select", y);
		selectWindow("C1-"+title1+"");
		run("Duplicate...", " ");
		setAutoThreshold("Shanbhag dark");
		setThreshold(992, 65535);
		run("Convert to Mask");
		run("Open");
		roiManager("Select", y);
		run("Analyze Particles...", "size=20-1000 circularity=0.10-1.00 show=Outlines exclude include add");	
	}
	print("saving QC image");
	selectWindow(""+title2+"-1.czi");
	Stack.setDisplayMode("composite");
	print("adjusting output B/C");
	run("Enhance Contrast", "saturated=0.35");
	run("Next Slice [>]");
	run("Enhance Contrast", "saturated=0.35");
	selectWindow(""+title2+"-1.czi");
	Stack.setDisplayMode("composite");
	roiManager("Select", newArray());
	roiManager("Show None");
	roiManager("Show All");
	run("Flatten");
	saveAs("Tiff", ""+dir2+"_"+title1+"_QC.tif");
	print("_");
	print("analysis:");
	run("Clear Results");
	ROIc = roiManager("count");
	while(ROIc>nROI){
		roiManager("Select", ROIc-1);
		roiManager("Delete");
		ROIc = ROIc-1;
	}
	nROI = parseInt(roiManager("count"));
	print("analysing "+actimage+"");
	selectWindow("ResOUT");
	IJ.renameResults("Results");
	for (z=0; z<nROI; z++) {
		n=z+1;
		print("analysing cell "+n+" of "+ROIc+" in "+actimage+"");
		roiManager("Select", z);
		run("Fit Ellipse");
		run("Set Measurements...", "area feret's redirect=None decimal=4");
		run("Measure");
		selectWindow("C1-"+title1+"");
		run("Duplicate...", " ");
		setAutoThreshold("Shanbhag dark");
		setThreshold(992, 65535);
		run("Convert to Mask");
		run("Open");
		roiManager("Select", z);
		//run("Set Measurements...", "area redirect=None decimal=4");
		run("Analyze Particles...", "size=20-1000 circularity=0.10-1.00 show=Outlines exclude include add");
		ROIc2=roiManager("count");
		nNUC=0;
		if(ROIc!=ROIc2){
			nNUC=parseInt(ROIc2-ROIc);
			for(j=ROIc2; j>(ROIc); j--){
				roiManager("Select", (j-1));
				roiManager("Delete");
			}
		}
			else{
				nNUC=0;
			}
		print(""+nNUC+" nuclei detected");
		setResult("nuclei", nResults-1, nNUC);
		setResult("filename", nResults-1, actimage);
		updateResults();
		ccount=nResults-1;
		}
		run("Close All");
	}
	print("saving results from "+N+" images ("+ccount+" cells) to "+dir2+"_"+day+".xls");
	print("_");
	selectWindow("Results");
	saveAs("txt", ""+dir2+"_"+day+"_"+month+".xls");
	print(""+N+" Images analysed");
	print("see output data in Destination Folder: "+dir2+"");
	print("_");
	getDateAndTime(year, month, week, day, hour, min, sec, msec);
	print("Finished analysis of "+N+" Images at: "+day+"/"+month+"/"+year+" :: "+hour+":"+min+":"+sec+"");
	selectWindow("Log");
	saveAs("Text", ""+dir2+"/log_analysis_"+day+"-"+month+"-"+year+"_"+hour+"h"+min+"min.txt");
	showMessage("Report",""+N+" Images analysed - see output data in Destination Folder: "+dir2+"");
	while (isOpen("ResOUT")) {
		selectWindow("ResOUT");
		run("Close");
    } 
    while (isOpen("Results")) {
		selectWindow("Results");
		run("Close");
    } 
	run("Close All");
	setBatchMode(false);
}
//JW_23.05