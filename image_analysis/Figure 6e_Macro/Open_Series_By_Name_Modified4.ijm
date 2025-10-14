//This macro will create a selection on channel 2 based on a threshold method, then measure in channel 4
//it will only work in images of the series that contain in the title the string "seriesName"

setOption("ExpandableArrays", true);

run("Set Measurements...", "area mean standard min centroid perimeter shape integrated display redirect=None decimal=5");

dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);


seriesName = "20";


run("Bio-Formats Macro Extensions");

for(f=0; f<files.length; f++) {
	projectname=files[f];
	file=dir+files[f];
	if(endsWith(file, ".lif")) {

	print(dir+file);
	Ext.setId(file);
	
	Ext.getSeriesCount(nSeries);
	seriesToOpen = newArray;
	sIdx = 0;
	for(i = 0; i < nSeries; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(name);
		// IJ.log(name);
		if( startsWith(name, seriesName) == true ) {
			seriesToOpen[sIdx++] = i + 1;
		}
	}
	
		Array.print(seriesToOpen);
		for(s = 0; s < seriesToOpen.length; s++) {
			run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
			imageName=getTitle();
			selectWindow(imageName);
			
			run("Split Channels");
			selectWindow("C1-" + imageName);
			run("Close");
			selectWindow("C3-" + imageName);
			run("Close");
			
			
			selectWindow("C2-" + imageName);
			run("8-bit");
			run("Manual Threshold...", "min=160 max=255");
			run("Create Selection");
			wait(500);
			run("Create Mask");
			
			selectWindow("Mask");
			run("Dilate");
			
			run("Create Selection");
			wait(500);
			roiManager("Add");
			selectWindow("ROI Manager");
			roiManager("Save", dir + "C2-" + imageName + ".zip");
			selectImage("C2-" + imageName);
			roiManager("Select", 0);
			run("Make Inverse");
			wait(500);
			run("Measure");
			
			wait(500);
			
			
			selectImage("C2-" + imageName);
			selectImage("C2-" + imageName);
			selectImage("C2-" + imageName);
			run("Close");
			selectWindow("Mask");
			selectWindow("Mask");
			selectWindow("Mask");
			run("Close");
			

			

		
			selectWindow("C4-" + imageName);
			run("8-bit");

			setAutoThreshold("Moments dark no-reset");
			run("Create Selection");
			wait(500);
			run("Create Mask");
			//wait(500);
			selectWindow("Mask");
			
			nrois=roiManager("Count");
			if (nrois == 0) {
				roiManager("Open", dir + "C2-" + imageName + ".zip");
			}
			selectWindow("Mask");
			roiManager("Select", 0); 
			run("Clear", "slice");

			selectWindow("Mask");
			
			selectWindow("ROI Manager");

			run("Close");
			
			selectWindow("Mask");
			getDimensions(width, height, channels, slices, frames);
			selectWindow("Mask");
			run("Select None");
			run("Analyze Particles...", "size=5-Infinity pixel add");
			wait(500);
			newImage("Untitled", "8-bit white", width, height, 1);
			roiManager("Deselect");
			selectWindow("Untitled");
			roiManager("Fill");
			selectWindow("ROI Manager");
			run("Close");
			run("Create Selection");
			wait(500);
			roiManager("Add");
			selectWindow("ROI Manager");
			roiManager("Save", dir + "C4-" + imageName + ".zip");

			
			selectWindow("C4-" + imageName);
			roiManager("Select", 0);

			roiManager("Measure");
			wait(500);
			
			selectWindow("C4-" + imageName);
			selectWindow("C4-" + imageName);
			selectWindow("C4-" + imageName);
			run("Select None");
			saveAs("tiff", dir + "C4-" + imageName);
			wait(500);
			selectWindow("C4-" + imageName + ".tif");
			run("Close");
			selectWindow("Mask");
			selectWindow("Mask");
			selectWindow("Mask");
			run("Close");
			
			roiManager("Deselect");
			selectWindow("ROI Manager");
			roiManager("Save", dir + "C4-" + imageName + ".zip");
			
			selectWindow("ROI Manager");
			run("Close");
			selectWindow("Untitled");
			run("Close");

		
		}
		
		
		
		selectWindow("Results");
		saveAs("results", dir + projectname+ "_CompleteResults.txt");
		
		nr=nResults();
		for (nn =0; nn<nr; nn++) {
			alabel=getResultString("Label", nn);
			if(startsWith(alabel, "C4-")) {

				IJ.deleteRows(nn, nn) ;
				nn=nn-1;
				nr=nr-1;
			}	
		}
		selectWindow("Results");
		saveAs("results", dir + projectname+ "_C2-Results.txt");	
		
		selectWindow("Results");
		run("Close");
		run("Results... ", "open=[" +dir + projectname + "_CompleteResults.txt]");
		nr=nResults();
		for (nn =0; nn<nr; nn++) {
			alabel=getResultString("Label", nn);
			if(startsWith(alabel, "C2-")) {

				IJ.deleteRows(nn, nn) ;
				nn=nn-1;
				nr=nr-1;
			}	
		}
		selectWindow("Results");
		saveAs("results", dir + projectname+ "_C4-Results.txt");	
		selectWindow("Results");
		run("Close");
		
	}
}