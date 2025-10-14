//This macro segments using CellPose and applies corrections based on area and circularity

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//parameters for small ld correction: adipocytes will be removed if they are bigger and darker than these limits
maxareabr=100;
minmeanbr=187;
//maximum area and min cir: all segmentations above this area or below this circ will be removed
maxarea=400;
mincirc=0.4;
//parameters for abhorrent shapes: segmentations that are bigger and less circular than these will be removed
maxareac=300;
mincircc=0.6;
//////////////////////////////////////////////////////////////////////////////////////////////////////////







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
			

			finalname=file+ " - " + filenaming + " " + iseriesname + "_allcrops0.zip";
			print(finalname);
			
			if(!File.exists(finalname)) {
				
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
				
				sizex=width/ncropsx;
				sizey=height/ncropsy;
				
				
				for ( i=0; i<ncropsx; i++) {
					for ( j=0; j<ncropsy; j++) {
						selectWindow(image8p);
						makeRectangle(0+i*sizex, 0+j*sizey, sizex, sizey);
						run("Duplicate...", " ");
						imagecropp=getTitle();
						selectWindow(imagecropp);
						
						selectWindow(image8);
						makeRectangle(0+i*sizex, 0+j*sizey, sizex, sizey);
						run("Duplicate...", " ");
						imagecrop=getTitle();
						selectWindow(imagecrop);
	
						//first segmentation of small LDs, cellpose, seems to ignore the big ones
						if(!File.exists(dir + image +"_" + i + "_"+ j + "_small0.zip")) {
							run("Cellpose Advanced", "diameter=5 cellproba_threshold=-6.0 flow_threshold=0.7 anisotropy=1.0 diam_threshold=100.0 model=cyto nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
							run("Label image to ROIs", "rm=[RoiManager[size=303, visible=true]]");
							
							nsrois=roiManager("Count");
							if (nsrois >0 ){
								roiManager("Deselect");
								roiManager("Save", dir+image+"_" + i + "_"+ j + "_small0.zip");
							}
						} else {
							print("small ones already segmented, we open them");
							roiManager("Open", dir+image+"_" + i + "_"+ j + "_small0.zip");
						}
						
						nsrois=roiManager("Count");
						if (nsrois >0 ){
							
							newImage("masksmall", "8-bit white", sizex, sizey, 1);
							wait(100);
							print("masksmall made");
							selectWindow("masksmall");
							roiManager("Show None");
							roiManager("Show All");
							setForegroundColor(0, 0, 0);
							setBackgroundColor(255, 255, 255);
							selectWindow("masksmall");
							roiManager("Show All");
							roiManager("Fill");
							selectWindow("ROI Manager");
							run("Close");
					
							selectWindow("masksmall");
							saveAs("Tiff", dir+ title +"_" + i + "_"+ j +  "_masksmall0.tif");
						} else {
							print("zero small adipos");
						}
						print("bb");
						//segmentation big
						if(!File.exists(dir + image +"_" + i + "_"+ j +  "_big0.zip")) {
							selectWindow(imagecrop);

							print("c");

							run("Duplicate...", "ignore");
							imaget=getTitle();
							selectWindow(imaget);
							setAutoThreshold("Default no-reset");

							run("Convert to Mask");

							wait(4000);
							run("Measure");
							wait(4000);
							rawint=getResult("RawIntDen", 0);
							area=getResult("Area", 0);
							fr=rawint/area;
							selectWindow("Results");
							wait(4000);
							run("Close");
							selectWindow(imaget);
							run("Close");
							
							selectWindow(imagecrop);
							run("Scale...", "x=0.5 y=0.5 interpolation=Bilinear average create");
							image025=getTitle();
							selectWindow(image025);
							
							if(fr>880) {
								percb=10;
							} else {
								percb=fr*0.0182-6.0435;
								if (percb<1) {
									percb=1;
								}
							}
							
							print("rawint: "+rawint);
							print("area: "+area);
							print("fr: "+fr);
							print("percb: "+percb);
							////////////
							run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'" + image025 + "', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'" + percb + "', 'percentileTop':'100.0', 'probThresh':'0.01', 'nmsThresh':'0.0', 'outputType':'ROI Manager', 'nTiles':'16', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'true', 'showCsbdeepProgress':'true', 'showProbAndDist':'false'], process=[false]");
							
							nbrois=roiManager("Count");
							if (nbrois >0 ){
								RoiManager.scale(2, 2, false);
								selectWindow(image);
								roiManager("Show None");
								roiManager("Show All");
								roiManager("Deselect");
								roiManager("Save", dir+image+"_" + i + "_"+ j + "_big0.zip");
							}
							selectWindow(image025);
							run("Close");
						}else {
								print("bigs already found");
								roiManager("Open", dir+image+"_" + i + "_"+ j + "_big0.zip");
								wait(1000);
								}
						

							print("d");
							nbrois=roiManager("Count");
							if (nbrois >0 ){
							if(!File.exists(dir + image +"_" + i + "_"+ j +  "_big1.zip")) {
								selectWindow(imagecropp);
								roiManager("Deselect");
								roiManager("Measure");
				
								nresults=nResults();
								print(nresults);
								
								errors=newArray(0);
								for (k=0; k<nResults; k++) {
									area=getResult("Area", k);
									mean=getResult("Mean", k);
									circ=getResult("Circ.", i);
									if( (area > maxareabr && mean < minmeanbr) || (area> maxarea) || (area > maxareac && circ < mincircc) || (circ < mincirc) ) { 
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
								
								roiManager("deselect");
								
								nbrois1=roiManager("Count");
								if (nbrois1 >0 ){
								roiManager("Save", dir+image+"_" + i + "_"+ j + "_big1.zip");
								}
							
							} else {
								print("bigs already corrected");
								roiManager("Open", dir+image+"_" + i + "_"+ j + "_big1.zip");
								wait(1000);
							}				
				
				
							
							
							print("dd");
							newImage("maskbig", "8-bit white", sizex, sizey, 1);
							wait(100);
							print("maskbig made");
							selectWindow("maskbig");
							nnrois1=roiManager("Count");
							if (nnrois1 >0 ){
							roiManager("Show None");
							roiManager("Show All");
							setForegroundColor(0, 0, 0);
							setBackgroundColor(255, 255, 255);
							selectWindow("maskbig");
							roiManager("Show All");
							roiManager("Fill");
							}

					
							selectWindow("maskbig");
							
							saveAs("Tiff", dir+ title +"_" + i + "_"+ j +  "_maskbig1.tif");
							titlemaskbig=getTitle();
						} else {
							print("zero big adipos");
						}
						
						if (nbrois>0 && nsrois > 0) {
						print("We merge now small and big segmentations");
					
						if(!File.exists(dir + image +"_" + i + "_"+ j +  "_all0.zip")) {
							print("Correcting small detections");

							selectWindow("ROI Manager");
							run("Close");
							selectWindow(title +"_" + i + "_"+ j +  "_maskbig1.tif");							
							roiManager("Open", dir+image+"_" + i + "_"+ j + "_small0.zip");
							selectWindow(titlemaskbig);
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
							roiManager("Save", dir+image+"_" + i + "_"+ j + "_smallcorrected0.zip");
							selectWindow(titlemaskbig);
							run("Close");

							if(File.exists(dir + image +"_" + i + "_"+ j +  "_big1.zip")) {
							roiManager("Open", dir+image +"_" + i + "_"+ j +  "_big1.zip");
							}
							roiManager("Deselect");
							RoiManager.translate(i*sizex, j*sizey);
							roiManager("Save", dir+image+"_" + i + "_"+ j + "_all0.zip");
							selectWindow("ROI Manager");
							run("Close");

						} else {
								print("this crop seems to be done aready");
								selectWindow("ROI Manager");
								run("Close");
							}
						selectWindow(imagecrop);
						run("Close");
						} else {
							print("no adipocytes found, image must be wrong");
						}
					}
				}
				
				//we have segmented all the crops, we now merge all segmentations in same image
				if (isOpen("ROI Manager")){
					selectWindow("ROI Manager");
					run("Close");
				}
								
				selectWindow(image8);
				
				for ( i=0; i<ncropsx; i++) {
					for ( j=0; j<ncropsy; j++) {
						if(File.exists(dir+image +"_" + i + "_"+ j +  "_all0.zip")){
						roiManager("Open", dir+image +"_" + i + "_"+ j +  "_all0.zip");
						}
					}
				}
				ntrois=roiManager("Count");
				if (ntrois>0){
				roiManager("Deselect");
				roiManager("Save", dir+image+ "_allcrops0.zip");
				}
				selectWindow("ROI Manager");
				run("Close");
				wait(500);
				run("Close All");
			
			}
			
		}
	}
}
		
