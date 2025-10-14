//This macro will use the user defined parameters to segment based on threshold in channels 2, 3 and 4, then measure in all the channels and calculate colocalization between segmentations
//it will also let you segment manually the regions containing eritrocytes, which are very bright, so they do not distrub the segmentation
//it can also open eritrocyte regions previously segmented manually
//it will also segment nuclei


////////////////////////////////////////////
//User defined parameters
seriesName = "all";
C3method="Triangle"; 
C4method="Moments"; 
C2method="Moments";
C3minsize="30"; 
C4minsize= "10";
C2minsize= "50";
C3dilates= 4;
C4dilates= 2;
C2dilates= 2;
minthreshold = 22;
maxthreshold = 255;
system="mac";
////////////////////////////////////////////


//setBackgroundColor(0, 0, 0);

setOption("ExpandableArrays", true);

run("Set Measurements...", "area mean standard min centroid perimeter shape integrated display redirect=None decimal=5");

dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);


run("Bio-Formats Macro Extensions");

//We are going to loop for each element of the series, do a max projection when necessary, segment channels 2-4 and nuclei in channel 1
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
		if(seriesName != "all") {

			if( indexOf(name, seriesName) != -1) {
				seriesToOpen[sIdx++] = i + 1;
			}
			} else { //we open all the series
				seriesToOpen[sIdx++] = i + 1;
			}
		}


	
	Array.show(seriesToOpen);
		

		for(s = 0; s < seriesToOpen.length; s++) {

			run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
			imageName=getTitle();
			selectWindow(imageName);
				
				if(!File.exists(dir +"C1-" +  imageName + "_Nuclei.zip") ){
					print(imageName);
					selectWindow(imageName);
					stacking(imageName, dir); //stacking is a function that will split channels if only one image open, and do z stack and split channels if several images open
					thresholding ("C3", C3method, C3minsize, C3dilates, imageName);
					thresholding ("C4", C4method, C4minsize, C4dilates, imageName);
					thresholding ("C2", C2method, C2minsize, C2dilates, imageName);
					print("Finding nuclei");
					selectWindow("C1" + "-" + imageName + ".tif");	
					
					run("Variance...", "radius=1");
                    run("Mean...", "radius=4");
                    run("Cellpose ...", "env_path=/Users/haiping/miniconda3/envs/cellpose env_type=conda model=nuclei model_path=/Applications/Fiji.app diameter=17 ch1=1 ch2=-1 additional_flags=[--cellprob_threshold, -1.0, --flow_threshold, 0.9]");
                    run("Label image to ROIs", "rm=[RoiManager[size=2782, visible=true]]");
					
					wait(2000);
					roiManager("Measure");
					selectWindow("Results");
					saveAs("results", dir + "C1-" + imageName + "_Nuclei.txt");
					roiManager("Save", dir +"C1-" +  imageName + "_Nuclei.zip");
					selectWindow("ROI Manager");
					run("Close");
					selectWindow("Results");
					run("Close");
				}
				//adipochasing is the function with all the instructions, I define it at the end
				run("Close All");
		} //end loop through all the series of each lif file		
	} // end condition if ends with lif
} //end loop through all the lifs in the folder



/////WE WILL ADD A CORRECTION TO TRY AND REMOVE ERITROCYTES
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
		if(seriesName != "all") {

			if( indexOf(name, seriesName) != -1) {
				seriesToOpen[sIdx++] = i + 1;
			}
			} else { //we open all the series
				seriesToOpen[sIdx++] = i + 1;
			}
		}



	
	Array.show(seriesToOpen);

		for(s = 0; s < seriesToOpen.length; s++) {
			run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
			imageName=getTitle();
			selectWindow(imageName);
			
			if(!File.exists(dir + "C4" + "-" + imageName + "corrected.zip") ){
			if(!File.exists( dir + imageName + "_Eritros.zip") ){
				getDimensions(width, height, channels, slices, frames);
				roiManager("Open", dir +"C1-" +  imageName + "_Nuclei.zip");
				newImage("Untitled", "8-bit white", width, height, 1);
				roiManager("Deselect");
				roiManager("Fill");
				selectImage("Untitled");
				run("Invert");

				
				selectWindow("ROI Manager");
				run("Close");
				
				roiManager("Open", dir + "C3" + "-" + imageName + ".zip");
				roiManager("Open", dir + "C2" + "-" + imageName + ".zip");
				roiManager("Open", dir + "C4" + "-" + imageName + ".zip");
				roiManager("Select", newArray(0,1,2));
				roiManager("AND");
				roiManager("Add");
				
				
				roiManager("Deselect");
				roiManager("Select", newArray(0,1,2));
				roiManager("Delete");
				
				nrois3=roiManager("Count");
				if(nrois3>0){
					roiManager("Select", 0);
					roiManager("Split");
					roiManager("Select", 0);
					roiManager("Delete");
					
					nrois=roiManager("Count");

					selectImage("Untitled");
					roiManager("Measure");
					
					dadipos=newArray(0);
					for ( ii=0; ii<nrois; ii++){
						rawintden=getResult("RawIntDen", ii);
						area=getResult("Area", ii);

						frcovered=(rawintden/255)/area;
						if(area < 150 || frcovered > 0.80){
							dadipos=Array.concat(dadipos, ii);
						}
					}
					roiManager("Deselect");

					
					selectWindow("Results");
					run("Close");
					
					roiManager("Deselect");
					selectWindow(imageName);
					Stack.setChannel(2);
					roiManager("Measure");
					for ( ii=0; ii<nrois; ii++){
						std=getResult("StdDev", ii);
						if(std < 10){
							dadipos=Array.concat(dadipos, ii);
						}
					}

					
					if (dadipos.length > 0) {
						roiManager("Select", dadipos);
						roiManager("Delete");
					}

					
					selectWindow("Results");
					run("Close");
					
					selectWindow("Untitled");
					run("Close");
					roiManager("Deselect");
					roiManager("Show All");
					
					waitForUser("Correct segmentations to exclude");
					roiManager("Deselect");

					nrois5=roiManager("Count");
					if(nrois5>0){
						roiManager("Deselect");
						roiManager("Combine");
						wait(1500);
						roiManager("Add");
						nrois2=roiManager("Count");
						nrois2=nrois2-1;
						roiManager("Select", nrois2);
						roiManager("save selected", dir + imageName + "_Eritros.zip");

						
						selectWindow("ROI Manager");
						run("Close");
						
						//Now we will remove this selection from the segmentations of channels 2, 3 and 4
						for (ch=2; ch<5; ch++) {
							roiManager("Open", dir + imageName + "_Eritros.zip");
							roiManager("Open", dir + "C" + ch + "-" + imageName + ".zip");
							roiManager("Select", newArray(0,1));
							roiManager("Combine");
							roiManager("Add");
							roiManager("Select", newArray(0,2));
							roiManager("XOR");
							roiManager("Add");
							roiManager("Select", 3);
							roiManager("save selected", dir + "C" + ch + "-" + imageName + "corrected.zip");
							selectWindow("ROI Manager");
							run("Close");
						}
					} else {
						for (ch=2; ch<5; ch++) {
							roiManager("Open", dir + "C" + ch + "-" + imageName + ".zip");
							roiManager("Save", dir + "C" + ch + "-" + imageName + "corrected.zip");
							selectWindow("ROI Manager");
							run("Close");
						}
					}
			
				} else {
					print("no overlaping between the three channes, no eritros found");
					waitForUser("No eritros found, press ok or draw eritrocyte areas");
					roiManager("Deselect");
					nrois4=roiManager("Count");
					if(nrois4>0){
						roiManager("Save", dir + imageName + "_Eritros.zip");
					} else {
						for (ch=2; ch<5; ch++) {
							roiManager("Open", dir + "C" + ch + "-" + imageName + ".zip");
							roiManager("Save", dir + "C" + ch + "-" + imageName + "corrected.zip");
							selectWindow("ROI Manager");
							run("Close");
						}
					}
				}
				run("Close All");
			} else {
				for (ch=2; ch<5; ch++) {
					roiManager("Open", dir + imageName + "_Eritros.zip");
					roiManager("Open", dir + "C" + ch + "-" + imageName + ".zip");
					roiManager("Select", newArray(0,1));
					roiManager("Combine");
					roiManager("Add");
					roiManager("Select", newArray(0,2));
					roiManager("XOR");
					roiManager("Add");
					roiManager("Select", 3);
					roiManager("save selected", dir + "C" + ch + "-" + imageName + "corrected.zip");
					selectWindow("ROI Manager");
					run("Close");
				}
				run("Close All");
			}
			
			
			}
			run("Close All"); 
		} //end loop through all the series of each lif file		
	} // end condition if ends with lif
} //end loop through all the lifs in the folder




//now we will loop through every tiff file of each channel, open the segmentations and measure in the other channels

ImageDir=dir;
list = getFileList(dir);

//define arrays
labels=newArray(0);

int3=newArray(0); int2=newArray(0); int4=newArray(0); int1=newArray(0);
area3=newArray(0); area2=newArray(0); area4=newArray(0);
seg3_int3=newArray(0); seg2_int2=newArray(0); seg4_int4=newArray(0);

seg3_int2=newArray(0); seg3_int4=newArray(0);
seg2_int3=newArray(0); seg2_int4=newArray(0);
seg4_int2=newArray(0); seg4_int3=newArray(0);

area3_coloc2=newArray(0); area3_coloc4=newArray(0);
area2_coloc3=newArray(0); area2_coloc4=newArray(0);
area4_coloc3=newArray(0); area4_coloc2=newArray(0);

area3_coloc24=newArray(0); area24_coloc3=newArray(0); 
		
		
for (n = 0; n < list.length; n++){	
	
	//CHANNEL 3
	if (endsWith(list[n], ".tif") == true && ( startsWith(list[n], "MAX_C3-") == true || startsWith(list[n], "C3-") == true)){
		open(dir+list[n]);
		imageName=getTitle();
		
		id=substring(imageName, indexOf(imageName, "C3-")+3); 
		//print(id);	
		prefix=substring(imageName, 0, indexOf(imageName, "C3-")+1); 
		print(prefix);
		labels=Array.concat(labels, id);
		id2=substring(id, 0, indexOf(id, ".tif")); 
		print(id2);
		
		run("Measure");
		selectWindow("Results");
		int_3=getResult("Mean", 0);
		selectWindow("Results");
		run("Close");
		
		roiManager("open", dir+ prefix +"3-"+id2+"corrected.zip");
		roiManager("open", dir+ prefix +"2-"+id2+"corrected.zip");
		roiManager("open", dir+ prefix +"4-"+id2+"corrected.zip");
		roiManager("deselect");
		roiManager("Measure");
		area_3=getResult("Area", 0);
		area_2=getResult("Area", 1);
		area_4=getResult("Area", 2);
		seg3__int3=getResult("Mean", 0);
		if(area_3==0) {
			seg3__int3=NaN;
			int_3=NaN;
		}
		seg2__int3=getResult("Mean", 1);
		if(area_2==0) {
			seg2__int3=NaN;
		}
		seg4__int3=getResult("Mean", 2);
		if(area_4==0) {
			seg4__int3=NaN;
		}
		int3=Array.concat(int3, int_3);
		
		area3=Array.concat(area3, area_3);
		area2=Array.concat(area2, area_2);
		area4=Array.concat(area4, area_4);
		seg3_int3=Array.concat(seg3_int3, seg3__int3);
		seg2_int3=Array.concat(seg2_int3, seg2__int3);
		seg4_int3=Array.concat(seg4_int3, seg4__int3);
		
		selectWindow("Results");
		run("Close");
		
		print("finding colocalization with FAP");
		nrprev=roiManager("count");
		roiManager("Select", newArray(0,1)); //coloc 3 and 2
		roiManager("AND");
		
		if(selectionType > -1) {
			roiManager("add");
		}
		nrpost=roiManager("count");
		if( nrpost>nrprev) {
			roiManager("select", nrpost-1);
			roiManager("Measure");
			coloc3_2=getResult("Area", 0);
			selectWindow("Results");
			run("Close");
		} else {
			coloc3_2=0;
		}
		
		
		print("finding colocalization of fib markers");
		nrprev=roiManager("count");
		roiManager("Select", newArray(0,2)); //coloc 3 and 4
		roiManager("AND");
		if(selectionType > -1) {
			roiManager("add");
		}
		nrpost=roiManager("count");
		if( nrpost>nrprev) {
			roiManager("select", nrpost-1);
			roiManager("Measure");
			coloc3_4=getResult("Area", 0);
			selectWindow("Results");
			run("Close");
		} else {
			coloc3_4=0;
		}
		

		print("finding colocalization with aSMA");
		nrprev=roiManager("count");
		roiManager("Select", newArray(1,2)); //coloc 2 and 4
		roiManager("AND");
		if(selectionType > -1) {
			roiManager("add");
		}
		nrpost=roiManager("count");
		if( nrpost>nrprev) {
			roiManager("select", nrpost-1);
			roiManager("Measure");
			coloc2_4=getResult("Area", 0);
			selectWindow("Results");
			run("Close");
		} else {
			coloc2_4=0;
		}
		
		print("finding colocalization with both fib markers");
		nrprev=roiManager("count");
		roiManager("Select", newArray(1,2)); //coloc 2 and 4
		roiManager("AND");
		if(selectionType > -1) {
			roiManager("add");
		}
		nrpost=roiManager("count");
		if( nrpost>nrprev) {
			
			roiManager("select", nrpost-1);
			roiManager("Measure");
			area_24=getResult("Area", 0);
			selectWindow("Results");
			run("Close");
			
			roiManager("Select", newArray(0,nrpost-1));
			roiManager("AND");
			if(selectionType > -1) {
			roiManager("add");
			}
			nrpost2=roiManager("count");
				
			if( nrpost2>nrpost) {
				roiManager("select", nrpost2-1);
				roiManager("Measure");
				coloc3_24=getResult("Area", 0);
				selectWindow("Results");
				run("Close");
			} else {
				coloc3_24=0;
				}
		} else {
			coloc3_24=0;
			area_24=0;
		}

		
		area3__coloc2 = coloc3_2 / area_3;
		area3__coloc4 = coloc3_4 / area_3;
		
		area2__coloc3 = coloc3_2 / area_2;
		area2__coloc4 = coloc2_4 / area_2;
		
		area4__coloc3 = coloc3_4 / area_4;
		area4__coloc2 = coloc2_4 / area_4;	
		
		area3__coloc24 = coloc3_24 / area_3;
		area24__coloc3 = coloc3_24 / area_24;
		
		area3_coloc2=Array.concat(area3_coloc2, area3__coloc2);
		area3_coloc4=Array.concat(area3_coloc4, area3__coloc4);
		area2_coloc3=Array.concat(area2_coloc3, area2__coloc3);
		area2_coloc4=Array.concat(area2_coloc4, area2__coloc4);
		area4_coloc3=Array.concat(area4_coloc3, area4__coloc3);
		area4_coloc2=Array.concat(area4_coloc2, area4__coloc2);
		
		area3_coloc24=Array.concat(area3_coloc24, area3__coloc24);
		area24_coloc3=Array.concat(area24_coloc3, area24__coloc3);
		
		run("Close All");
		
		
		open(dir+ prefix + "2-" + id);
		imageName2=getTitle();
		run("Measure");
		selectWindow("Results");
		int_2=getResult("Mean", 0);
		selectWindow("Results");
		run("Close");
		
		roiManager("deselect");
		roiManager("Measure");
		seg3__int2=getResult("Mean", 0);
		seg2__int2=getResult("Mean", 1);
		seg4__int2=getResult("Mean", 2);
		
		
		if(area_3==0) {
			seg3__int2=NaN;
		}
		if(area_2==0) {
			seg2__int2=NaN;
			int_2=NaN;
		}
		if(area_4==0) {
			seg4__int2=NaN;
		}
		int2=Array.concat(int2, int_2);
		
		seg3_int2=Array.concat(seg3_int2, seg3__int2);
		seg2_int2=Array.concat(seg2_int2, seg2__int2);
		seg4_int2=Array.concat(seg4_int2, seg4__int2);
		
		selectWindow("Results");
		run("Close");
		run("Close All");
		
		
		
		open(dir+ prefix + "4-" + id);
		imageName4=getTitle();
		run("Measure");
		selectWindow("Results");
		int_4=getResult("Mean", 0);
		selectWindow("Results");
		run("Close");
		
		roiManager("deselect");
		roiManager("Measure");
		seg3__int4=getResult("Mean", 0);
		seg2__int4=getResult("Mean", 1);
		seg4__int4=getResult("Mean", 2);
		
		if(area_3==0) {
			seg3__int4=NaN;
		}
		if(area_2==0) {
			seg2__int4=NaN;
		}
		if(area_4==0) {
			seg4__int4=NaN;
			int_4=NaN;
		}
		
		int4=Array.concat(int4, int_4);
		
		seg3_int4=Array.concat(seg3_int4, seg3__int4);
		seg2_int4=Array.concat(seg2_int4, seg2__int4);
		seg4_int4=Array.concat(seg4_int4, seg4__int4);
		
		selectWindow("Results");
		run("Close");
		run("Close All");
		
		
		selectWindow("ROI Manager");
		run("Close");
		
		open(dir+ prefix + "1-" + id);
		imageName1=getTitle();
		run("Measure");
		selectWindow("Results");
		int_1=getResult("Mean", 0);
		int1=Array.concat(int1, int_1);
		selectWindow("Results");
		run("Close");
		
		run("Close All");
		
	}
	
Array.show(labels, int3, int2, int4, int1, area3, area2,area4, 
	seg3_int3, seg2_int2, seg4_int4,
	seg3_int2, seg3_int4, seg2_int3, seg2_int4, seg4_int2, seg4_int3,
	area3_coloc2, area3_coloc4, area2_coloc3, area2_coloc4, area4_coloc3, area4_coloc2,
	area3_coloc24, area24_coloc3);
	
	
}

Array.show(labels, int3, int2, int4, int1, area3, area2,area4, 
	seg3_int3, seg2_int2, seg4_int4,
	seg3_int2, seg3_int4, seg2_int3, seg2_int4, seg4_int2, seg4_int3,
	area3_coloc2, area3_coloc4, area2_coloc3, area2_coloc4, area4_coloc3, area4_coloc2,
	area3_coloc24, area24_coloc3);

selectWindow("Arrays");

saveAs("results", dir +  "_FinalResults.txt");

		





































function stacking(nameofimage, directory) {
	
	titleList = getList("image.titles"); 
	titlecount = lengthOf(titleList);
	    
	if (titlecount > 1) {
		run("Concatenate...", "all_open open");
		selectWindow("Untitled");
		rename(nameofimage);
		getDimensions(width, height, channels, slices, frames);
		run("Split Channels");
		channels=channels+1;
		for (h=1; h<channels; h++) {
			selectWindow("C" + h + "-"+nameofimage);
			run("Z Project...", "projection=[Max Intensity]");
			selectWindow("C" + h + "-"+nameofimage);
			run("Close");
		}
		imgs = getList("image.titles");
		for (i = 0; i < imgs.length; i++) {
			selectImage(imgs[i]);
   			theimagename=getTitle();
   			saveAs("tiff", directory + theimagename);
		}
		
	} else {
		
		if (titlecount > 0) {
			print("no images open");
			getDimensions(width, height, channels, slices, frames);
			run("Split Channels");
			imgs = getList("image.titles");
			for (i = 0; i < imgs.length; i++) {
				selectImage(imgs[i]);
	   			theimagename=getTitle();
	   			saveAs("tiff", directory + theimagename);
			}
		}
	}
}



function thresholding(channel, method, size, ndilates, imagename) {
	selectWindow(channel + "-" + imagename + ".tif");	
	run("8-bit");
	
	if(method != "manual") {
		setAutoThreshold(method + " dark no-reset");
	} else {
		setThreshold(minthreshold, maxthreshold, "raw");
	}

	getDimensions(width, height, channels, slices, frames);
	run("Create Mask");

	run("Analyze Particles...", "size=" + size +"-Infinity pixel add");
	roiManager("Show All without labels");
	wait(200);
	print("particles analysed");
	nroisi=roiManager("count");
	if (nroisi > 0) {
		roiManager("deselect");
		roiManager("Combine");
		roiManager("Add");
		nroism=roiManager("count");
		nroism=nroism-1;
		roiManager("select", nroism);
		if (system== "mac") {
		setBackgroundColor(255, 255, 255);
		} else {
			setBackgroundColor(0, 0, 0);
		}
		run("Clear Outside");
			

		roiManager("Show None");
		roiManager("Show All");

		wait(2000);
		selectWindow("ROI Manager");
		wait(500);
		run("Close");
		selectWindow("mask");
		print("doing dilates");
		
		counter=0;
		while (counter < ndilates) {
			selectWindow("mask");
			run("Dilate");
			counter=counter+1;
		}
		run("Create Selection");
		wait(2000);	
		
		roiManager("Add");
		wait(2000);
		print("selection added to roi manager");
		nc=roiManager("Count");
		print(nc);
		while(nc<1){
			roiManager("Add");
			nc=roiManager("Count");
		}
		wait(2000);
		
		print("selection amplified");
	
	
	} else { 
		print("no positive area detected");
		makePoint(0,0);
		roiManager("Add");
		wait(2000);
	}
	

	
	selectWindow("mask");
	selectWindow("mask");
	selectWindow("mask");
	selectWindow("mask");
	run("Close");

	selectWindow("ROI Manager");
	wait(200);
	roiManager("Save", dir + channel + "-" + imagename + ".zip");
	print("selection saved");
	selectWindow("ROI Manager");
	run("Close");
	//}

	selectWindow(channel + "-" + imagename+ ".tif");
	resetThreshold();
}






