//This macro segments specifically LDs of very small area, indicated in tissues showing browning

getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("started " + hour + ":" + minute + ":" + second);

setOption("ExpandableArrays", true);

run("Set Measurements...", "area mean standard min centroid perimeter shape integrated display redirect=None decimal=5");

dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);

resolution=3;
seriesName = "all";
ncropsx= 2;
ncropsy=2;


run("Bio-Formats Macro Extensions");

for(f=0; f<files.length; f++) {
	projectname=files[f];
	file=dir+files[f];
	if(endsWith(file, ".czi")) {

	print(dir+file);
	Ext.setId(file);
	
	Ext.getSeriesCount(nSeries);
	seriesToOpen = newArray;
	sIdx = 0;
	sizeXprev=0;
	for(i = 0; i < nSeries; i++) {
		Ext.setSeries(i);
		Ext.getSeriesName(name);
		Ext.getSeriesMetadataValue("SizeX", value);
		Ext.getSizeX(sizeX);
		// IJ.log(name);
			if(seriesName != "all") {

				if( indexOf(name, seriesName) != -1) {
					if (sizeX > sizeXprev) {
					ii=i+(resolution-1);
					if (ii < nSeries) {
						seriesToOpen[sIdx++] = ii + 1;
					}

					}
					sizeXprev=sizeX;
				}
			} else { //we open all the series
				if (sizeX > sizeXprev) {
				ii=i+(resolution-1);
				if (ii < nSeries) {
					seriesToOpen[sIdx++] = ii + 1;
				}

				}
				sizeXprev=sizeX;
			}		
	}
	
		Array.print(seriesToOpen);
		for(s = 0; s < seriesToOpen.length; s++) {
			print(file);
			where=lastIndexOf(file, "\\");
			where=where+1;
			filenaming=substring(file, where);

			if (seriesToOpen[s] > 9) {
				iseriesname="#"+seriesToOpen[s];
			}else {
				iseriesname="#0"+seriesToOpen[s];
			}
			print("iseriesname: "+ iseriesname);
			

			finalname=file+ " - " + filenaming + " " + iseriesname + "__allcrops1.zip";
			print(finalname);
			
			if(!File.exists(finalname)) {
				//create mask of the already found LDs
				
				run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
				imageName=getTitle();
				image=imageName;
				title=image;
				selectWindow(image);
				
				run("RGB Color");

				run("8-bit");
				image8p=getTitle();
				run("Duplicate...", " ");

				run("Enhance Contrast...", "saturated=4 equalize"); //NEW in 5

				image8=getTitle();
				selectWindow(image8);
				getDimensions(width, height, channels, slices, frames);
				selectWindow(image8);
				
				newImage("maskall", "8-bit white", width, height, 1);
				roiManager("Open", file+ " - " + filenaming + " " + iseriesname + "_allcrops0.zip");
				
				
				
				selectWindow("maskall");
				roiManager("Show None");
				roiManager("Show All");
				setForegroundColor(0, 0, 0);
				setBackgroundColor(255, 255, 255);
				selectWindow("maskall");
				roiManager("Show All");
				roiManager("Fill");
				selectWindow("ROI Manager");
				run("Close");
				
				sizex=width/ncropsx;
				sizey=height/ncropsy;
				
				
				for ( i=0; i<ncropsx; i++) {
					for ( j=0; j<ncropsy; j++) {
						selectWindow(image8);
						makeRectangle(0+i*sizex, 0+j*sizey, sizex, sizey);
						run("Duplicate...", " ");
						imagecrop=getTitle();
						selectWindow(imagecrop);
	
						//first segmentation of small LDs, cellpose, seems to ignore the big ones
						if(!File.exists(dir + image +"_" + i + "_"+ j + "_tiny0.zip")) {
							run("Cellpose Advanced", "diameter=3 cellproba_threshold=-6.0 flow_threshold=2 anisotropy=1.0 diam_threshold=200.0 model=cyto nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
							run("Label image to ROIs", "rm=[RoiManager[size=303, visible=true]]");
							
							nsrois=roiManager("Count");
							if (nsrois >0 ){
								roiManager("Deselect");
								RoiManager.translate(i*sizex, j*sizey);
								roiManager("Save", dir+image+"_" + i + "_"+ j + "_tiny0.zip");
							}
						} else {
							print("tiny ones already segmented, we open them");

						}
						
						selectWindow(imagecrop);
						run("Close");
						
						if(isOpen("ROI Manager")) {
							selectWindow("ROI Manager");
							run("Close");
						}
						
					}
				}
				
				for ( i=0; i<ncropsx; i++) {
					for ( j=0; j<ncropsy; j++) {
						if(File.exists(dir+image +"_" + i + "_"+ j +  "_tiny0.zip")){
						roiManager("Open", dir+image +"_" + i + "_"+ j +  "_tiny0.zip");
						}
					}
				}
				
				nsrois=roiManager("Count");
				if (nsrois >0 ){
					selectWindow("maskall");
					run("Invert");
					roiManager("Deselect");
					roiManager("Measure");
					nresults=nResults();
					print(nresults);
					
					errors=newArray(0);
					for (k=0; k<nResults; k++) {
						rawintden=getResult("RawIntDen", k);
						area=getResult("Area", k);
						propoverlap=(rawintden/255)/area;
						mean=getResult("Mean", k);
						if( propoverlap > 0.1) { //if they overlap
							errors= Array.concat(errors, k);
						}
					}
					selectWindow("Results");
					run("Close");

					//We remove big selections overlapping with the small ones				
					if( errors.length > 0) {
						roiManager("select", errors);
						roiManager("delete");
						roiManager("deselect");
					}
					
					selectWindow("maskall");
					run("Close");
					
					
					selectWindow(image8p);
					roiManager("Deselect");
					roiManager("Measure");
					nresults=nResults();
					print(nresults);
					
					errors=newArray(0);
					for (k=0; k<nResults; k++) {
						mean=getResult("Mean", k);

						if( mean > 180) { //if they overlap
							errors= Array.concat(errors, k);
						}
					}
					selectWindow("Results");
					run("Close");
					//We remove big selections overlapping with the small ones				
					if( errors.length > 0) {
						roiManager("select", errors);
						roiManager("delete");
						roiManager("deselect");
					}
					roiManager("Save",file+ " - " + filenaming + " " + iseriesname +  "_tinycorrected0.zip");
					
					
				} else {
					print("zero tiny adipos");
				}
				print("bb");
						
						
				//open segmentation big and small
				roiManager("Deselect");
				roiManager("Open", file+ " - " + filenaming + " " + iseriesname + "_allcrops0.zip");		
				roiManager("Save", dir+image+"_" + "_allcrops1.zip");		

				selectWindow("ROI Manager");
				run("Close");
				
				run("Close All");
			
			}
		}
	}
}
		
