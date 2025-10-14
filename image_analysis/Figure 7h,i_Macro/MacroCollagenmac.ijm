if (isOpen("Results")) {
	selectWindow("Results");
	run("Close");
}

run("Set Measurements...", "area mean standard min shape integrated limit display redirect=None decimal=5");
dir=getDirectory("Select folder with images");
list = getFileList(dir);

for (n=0; n<list.length; n++){
	if(endsWith(list[n], ".tiff") == true){
		
		open(dir+list[n]);
	
		title=getTitle();
		print(title);
		
		run("Colour Deconvolution", "vectors=[Masson Trichrome]");
		
		selectImage(title+"-(Colour_3)");
		close();
		selectImage("Colour Deconvolution");
		close();
		
		selectImage(title+"-(Colour_1)");
		getDimensions(width, height, channels, slices, frames);
		setAutoThreshold("Default no-reset");

		selectImage(title+"-(Colour_1)");
		check=is("binary");
		while(check<1){
			run("Convert to Mask");
			check=is("binary");
		}

		wait(250);
		print("Blue mask made");
		
		selectImage(title+"-(Colour_2)");

		setAutoThreshold("Default no-reset");
		run("Convert to Mask");
		run("Dilate");
		run("Create Selection");
		wait(250);
		
		if (selectionType() > -1){
			roiManager("Add");
			print("Red selection made");
			selectImage(title+"-(Colour_1)");
			roiManager("Select", 0);
			roiManager("Deselect");
			run("Clear", "slice");
		}else {
			print("no red thresolded signal");
			selectImage(title+"-(Colour_1)");
		}
		
		

		
		
		if (isOpen("ROI Manager")) {
			selectWindow("ROI Manager");
			run("Close");
		}
		
		run("Select None");
		
		
		print("analyse particles");
		selectImage(title+"-(Colour_1)");
		run("Open");
		run("Analyze Particles...", "size=5-Infinity exclude add");

		
		newImage("MaskParticles", "8-bit white", width, height, 1);
		nrois=roiManager("Count");
		print("particles analysed");
		
		if (nrois>0){
			roiManager("Deselect");
			roiManager("Fill");
		} else {
			print("no particles found in blue channel, no collagen");
			
		}
		
		print("particles filled");

		run("Select None");
		
		if (isOpen("ROI Manager")) {
			selectWindow("ROI Manager");
			run("Close");
		}
		
		selectImage("MaskParticles");
		run("Create Selection");
		wait(500);

		
		countrois=roiManager("Count");
		while(countrois<1){
			roiManager("Add");
			countrois=roiManager("Count");
			
		}
		
		wait(250);
		print("nrois="+countrois);

		
		
		nrois2=roiManager("Count");
		print("nrois="+nrois2);
		print("saving particles");
		
		selectWindow(title);
		if (nrois2>0){
			roiManager("Measure");
			
			roiManager("Save", dir+title+"_collagen.zip");
			selectWindow("ROI Manager");
			run("Close");
		} else {
			waitForUser("");
			print("no particles found in blue channel, me make a fake point to measure zero area");
			makePoint(0,0);
			roiManager("Add");
			roiManager("Measure");
			selectWindow("ROI Manager");
			run("Close");
		}
		
		

		run("Close All");		
	}

}

selectWindow("Results");
saveAs("Results", dir+"_ResultsCollagen.txt");

selectWindow("Results");
run("Close");

