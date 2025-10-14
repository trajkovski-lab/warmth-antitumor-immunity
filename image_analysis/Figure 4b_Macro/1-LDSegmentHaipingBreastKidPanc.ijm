////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is a macro to segment LDs, but it actually segments any round objects.
//It can use cellpose or stardist to segment
//and it has several options optimized for histological images, fluorescence and brightfield
//user must define parameters depending on the type of image and the expected diameter
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
extension=".png";
seriesfile=false;
methodbig="cellpose"; //either this or stardist
imagetype="histo"; //either fluo or histo. For brightfield images, put "fluo" too.
BF=false; //put true here if brightfield
channelld=1; //IMPORTANT. Here is the channel of whatever you want to segment
dosmall=true; // if you only want to segment cells, not LDs, or all your LDs are big (no browning expected) put false here and you save A LOT of time.
enhancecontrastfluo=false; // put true if you want the contrast of the fluo image to be increased. It will saturate a percentage of the pixels.
percsaturated=0.35; // percentage of pixels saturated if enhancecontrastfluo is set to true
blurimagefluo=false; //recommended false. Set to true if image too noisy, or too many small segmentations found. It sacrifices small LD segmentations, but sometimes it compensates
radiusblur=1; //the more the value, the more it blurs. Only works if blurimagefluor is true

//settings for cellpose segmentation of small and big lds
//optimaze for each experiment
diamsmall=70;
flowsmall=0.4;
probsmall=2.0;
modelsmall="cyto2"; //choose one of available cellpose models

diambig=200; 
flowbig=0.4;
probbig=2.0;
modelbig="cyto2"; 

//settings for find maxima to eliminate big LDs segmentations that correspond to merged small LDs
testmaxima=false;
topmaximapoints=2; // big selections with more maxima points than these, will be eliminated
maximaprominence=10;  //maxima detection sensitivity, smaller prominence will give more maxima points, increase if too many big LDs eliminated
dosigmamaxima=false; //should image be blurred before running maxima?
sigmamaxima=2; // by how much should it be blurred? only works if dosigmamaxima is true

//conditions for decisions between big and small ld segmentations
maxoverlapbigsmall=0.25; //how much big lds are allowed to overlap with small without being excluded
//parameters for small ld correction: adipocytes will be removed if they are bigger and darker than these limits
//area, intensity, size filters
maxareabr=999999999999999999999;
minmeanbr=0; //150;

//parameters to exclude big
//minimum intensity allowed
minmeanbig=190; //175 histo

//maximum area and min cir: all segmentations above this area or below this circ will be removed
maxarea=50000000000000000000000000000000000000000000000000000000000000000; 
mincirc=0.3;
//parameters for abhorrent shapes: segmentations that are bigger and less circular than these will be removed
maxareac=9999999999999999999999999999999999999999999999; 
minarea=1000;
maxmean=255;
minmean=190;

//parameters for image series and image cropping
resolution=3; //maginification you want to analyse. Higly recommend 3 (3rd best); If more resolution (1 or 2), the image will need to be cropped. Less resolution (>3) will greatly compromise small LD segmentation. If only big LDs expected, no problem with using >3
seriesName = "all"; //if set to "all", it does all the elemets of the series with the specified resolution. If you have a string for some set of images to analise, you can run macro on only those and exclude the rest
ncropsx= 1; //number of crops in x
ncropsy=1; //number of crops in y
if(ncropsx<1) {ncropsx=1;} //for caution
if(ncropsy<1) {ncropsy=1;} //for caution

measure=true; //true if we want to measure in the end. Put false if we are going to run the correction macro after
//////////////////////////////////////////////////////////////////////////////////////////////////////////







getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print("started " + hour + ":" + minute + ":" + second);

setOption("ExpandableArrays", true);
run("Set Measurements...", "area mean standard min centroid perimeter shape integrated display redirect=None decimal=5");

dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);

run("Bio-Formats Macro Extensions");

for(f=0; f<files.length; f++) {
	projectname=files[f];
	file=dir+files[f];
	if(endsWith(file, extension)) {

	print(dir+file);
	
	if (seriesfile == true) {
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
			
			ntotalfiles=seriesToOpen.length;
	} else {
		ntotalfiles = 1;
		}
	
	
	
		
		
		
		for(s = 0; s < ntotalfiles; s++) {
			print(file);
			
			if (seriesfile == true) {
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
				
				if(measure==true){
					finalname=file+ " - " + filenaming + " " + iseriesname + "_RESULTSFINAL.txt";
				}
			} else {
				finalname=file+  "_ch" +channelld+ "_allcrops0.zip";
				if(measure==true){
					finalname=file+  "_ch" +channelld+ "_RESULTSFINAL.txt";;
				}
			}
			
			print(finalname);
			
			if(!File.exists(finalname)) {
				
				if (seriesfile == true) {
					run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
				} else{
					if( BF==true){
						print("doing brightfield image, importing grayscale");
						open(file);
					} else {
						run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
					}
				}
				imageName=getTitle();
				image=imageName;
				title=image;
				selectWindow(image);
				
				if(imagetype == "histo") {
				run("RGB Color");

				run("8-bit");
				image8p=getTitle();
				run("Duplicate...", " ");

				run("Enhance Contrast...", "saturated=4 equalize"); //NEW in 5

				} else {
					if (imagetype =="fluo") {
						selectWindow(image);
						Stack.setChannel(channelld);
						run("Duplicate...", "duplicate channels=" + channelld);
						
						if( blurimagefluo==true){
							run("Mean...", "radius="+ radiusblur);
						}
						
						if( enhancecontrastfluo==true){
							run("Enhance Contrast...", "saturated=" + percsaturated + " equalize");
						}
						
						image8p=getTitle();
					} else {
						print("IMAGE TYPE MUST BE EITHER fluo OR histo");
					}
				}
				
				
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
						if (dosmall== true) {
							if(!File.exists(dir + image +"_" + i + "_"+ j +  "_ch" +channelld+"_small0.zip")) {
								
								if(BF==true){
									run("Console");
									run("Cellpose Advanced", "diameter=" + diamsmall +" cellproba_threshold=" + probsmall + " flow_threshold=" + flowsmall + " anisotropy=1.0 diam_threshold=100.0 model=" + modelsmall + " nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
									run("Label image to ROIs", "rm=[RoiManager[size=303, visible=true]]");
								} else {
									run("Console"); 
									run("Cellpose Advanced", "diameter=" + diamsmall +" cellproba_threshold=" + probsmall + " flow_threshold=" + flowsmall + " anisotropy=1.0 diam_threshold=100.0 model=" + modelsmall + " nuclei_channel=0 cyto_channel=0 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
									run("Label image to ROIs", "rm=[RoiManager[size=303, visible=true]]");									
								}

								
								nsrois=roiManager("Count");
								if (nsrois >0 ){
									roiManager("Deselect");
									roiManager("Save", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_small0.zip");
								}
							} else {
								print("small ones already segmented, we open them");
								roiManager("Open", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_small0.zip");
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
								saveAs("Tiff", dir+ title +"_" + i + "_"+ j +   "_ch" +channelld+ "_masksmall0.tif");
							} else {
								print("zero small adipos");
							}
						} else {
							nsrois=0;
							print("small objects ignored as per user settings");
						}
						
						
						print("bb");
						//segmentation big
						if(!File.exists(dir + image +"_" + i + "_"+ j +  "_ch" +channelld+ "_big0.zip")) {
							selectWindow(imagecrop);
							//
							print("c");
							if (methodbig=="stardist") {
								print("chosen method for big LDs: stardist");

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
							}

							
							
							selectWindow(imagecrop);
							if (methodbig=="stardist") {
								run("Scale...", "x=0.5 y=0.5 interpolation=Bilinear average create");
							}
							image025=getTitle();
							selectWindow(image025);
							
							if (methodbig=="stardist") {
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
							} else {
								if(methodbig=="cellpose") {
									print("chosen method for big LDs: cellpose");

									run("Cellpose Advanced", "diameter=" + diambig + " cellproba_threshold=" + probbig + " flow_threshold=" + flowbig +" anisotropy=1.0 diam_threshold=100.0 model=" + modelbig + " nuclei_channel=2 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
									run("Label image to ROIs", "rm=[RoiManager[size=303, visible=true]]");

								} else {
									print("NOTICE: methodbig must be eihter stardist or cellpose");
								}
							}


							nbrois=roiManager("Count");
							if (nbrois >0 ){
								if (methodbig=="stardist") {
									RoiManager.scale(2, 2, false);
								}
								
								selectWindow(image);
								roiManager("Show None");
								roiManager("Show All");
								roiManager("Deselect");
								roiManager("Save", dir+image+"_" + i + "_"+ j + "_ch" +channelld+ "_big0.zip");
							}
							selectWindow(image025);
							run("Close");
						}else {
								print("bigs already found");
								roiManager("Open", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_big0.zip");
								wait(1000);
								}
						

							print("Eliminating big selections depending on (maxima), area and intensity parameters");
							
							if(isOpen("Results")){
								selectWindow("Results");
								run("Close");
							}
							
							
							//big selections are open in the ROI Manager right now
							nbrois=roiManager("Count");
							if (nbrois >0 ){
								if(!File.exists(dir + image +"_" + i + "_"+ j +   "_ch" +channelld+"_big1.zip")) {
									selectWindow(imagecropp);
									roiManager("Deselect");
									roiManager("Measure");
					
									nresults=nResults();
									print(nresults);
									
									errors=newArray(0);
									for (k=0; k<nResults; k++) {
										area=getResult("Area", k);
										mean=getResult("Mean", k);
										circ=getResult("Circ.", k);
										if( (area > maxareabr && mean < minmeanbr) || (area< minarea) ||(area> maxarea) || (area > maxareac && circ < mincircc) || (circ < mincirc)  || (mean < minmeanbig ) || area> maxarea || area < minarea || mean>maxmean || mean>maxmean){ 
											errors= Array.concat(errors, k);
										}
									}

									selectWindow("Results");
									run("Close");
									
									if(testmaxima==true){
										roiManager("Deselect");
										print("test maxima parameter set to true");
										selectWindow(imagecropp);
										if(imagetype=="fluo") {
											run("Duplicate...", " ");
											imageblur=getTitle();
											if (dosigmamaxima == true) {
												run("Gaussian Blur...", "sigma=" + sigmamaxima);
											}
										}
										run("Find Maxima...", "prominence=" + maximaprominence +" output=[Single Points]");
										if(imagetype=="fluo") {
											selectWindow(imageblur);
											run("Close");
										}
										imagemaxima=getTitle();
										selectWindow(imagemaxima);

										roiManager("Deselect");
										//this is to delete maxima too close to the ld segmentation border, we thicken rois and paint them black on the top of the maxima mask
										roiManager("Set Line Width", 3);
										setForegroundColor(0, 0, 0);
										setBackgroundColor(255, 255, 255);

										roiManager("Draw");
										roiManager("Set Line Width", 0);
										
										roiManager("Deselect");
										roiManager("Measure");
										
										errorsmax=newArray(0);
										for (mm=0; mm<nResults; mm++) {
											mrawintden=getResult("RawIntDen", mm);
											if(mrawintden> topmaximapoints*255) {
												errorsmax= Array.concat(errorsmax, mm);
											}
										}

										roiManager("deselect");
										
										if(lengthOf(errors)>0){
											roiManager("Select", errors);
											roiManager("Set Color", "pink");
											roiManager("Set Line Width", 10);
										}
										
										if(lengthOf(errorsmax)>0){
											roiManager("Select", errorsmax);
											roiManager("Set Color", "red");
											roiManager("Set Line Width", 10);
										}
										selectWindow(imagemaxima);
										run("Close");

										selectWindow("Results");
										run("Close");
										
										errors= Array.concat(errors, errorsmax);
									}
																
									//We remove big selections that meet the error conditions			
									if( errors.length > 0) {
										roiManager("select", errors);
										roiManager("delete");
										roiManager("deselect");
									}
									
									roiManager("deselect");
									
									nbrois1=roiManager("Count");
									if (nbrois1 >0 ){
										
									if (ncropsx==1 && ncropsy==1 && dosmall==false){
										roiManager("Save", dir+image+  "_ch" +channelld+"_allcrops0.zip");
									} else {
									roiManager("Save", dir+image+"_" + i + "_"+ j + "_ch" +channelld+ "_big1.zip");
									}
									}
								
								} else {
									print("bigs already corrected");
									roiManager("Open", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_big1.zip");
									wait(1000);
								}				
					
					
								
								if (ncropsx!=1 || ncropsy!=1 || dosmall!=false){
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
									
									saveAs("Tiff", dir+ title +"_" + i + "_"+ j +   "_ch" +channelld+ "_maskbig1.tif");
									titlemaskbig=getTitle();
								}
						} else {
							print("zero big adipos");
						}
						
						if (ncropsx!=1 || ncropsy!=1 || dosmall!=false){ //what if there are several crops but no small, test that: PEDING
							if (nbrois>0 && nsrois > 0) {
								print("We merge now small and big segmentations");
							
								if(!File.exists(dir + image +"_" + i + "_"+ j +  "_ch" +channelld+ "_all0.zip")) {
									print("Correcting small detections");

									selectWindow("ROI Manager");
									run("Close");
									selectWindow(title +"_" + i + "_"+ j +  "_ch" +channelld+ "_maskbig1.tif");							
									roiManager("Open", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_small0.zip");
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
										if( propoverlap > maxoverlapbigsmall ) { //if they overlap
											errors= Array.concat(errors, k);
										}
									}
									selectWindow("Results");
									run("Close");
									
								
									
									//NEWWWW
									selectWindow(image8p);
									roiManager("Deselect");
									roiManager("Measure");
									nresults=nResults();
									print(nresults);
									for (k=0; k<nResults; k++) {
										area=getResult("Area", k);
										mean=getResult("Mean", k);
										if( area> maxarea || area < minarea || mean>maxmean || mean<minmean ) { //if they overlap
											errors= Array.concat(errors, k);
										}
									}
									selectWindow("Results");
									run("Close");
									
									//////////////////////////////////////
									
									print("check1");

									//We remove big selections overlapping with the small ones				
									if( errors.length > 0) {
										roiManager("select", errors);
										roiManager("delete");
										roiManager("deselect");
									}

									roiManager("Save", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_smallcorrected0.zip");
									selectWindow(titlemaskbig);
									run("Close");
	
									if(File.exists(dir + image +"_" + i + "_"+ j +  "_ch" +channelld+ "_big1.zip")) {
										roiManager("Open", dir+image +"_" + i + "_"+ j +  "_ch" +channelld+ "_big1.zip");
									}
									roiManager("Deselect");
									RoiManager.translate(i*sizex, j*sizey);
									roiManager("Save", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_all0.zip");
									selectWindow("ROI Manager");
									run("Close");

								} else {
										print("this crop seems to be done aready");
										selectWindow("ROI Manager");
										run("Close");
									}
									
								print("A");
								if(isOpen(imagecrop)){
									selectWindow(imagecrop);
									run("Close");
								}
		
								print("B");
							} else {
								if (nbrois>0 ) {
									if(File.exists(dir + image +"_" + i + "_"+ j +  "_ch" +channelld+ "_big1.zip")) {
										roiManager("Open", dir+image +"_" + i + "_"+ j +  "_ch" +channelld+ "_big1.zip");
									} else {
										print("this shouldn't happen, check macro, there is some error");
										error
									}
									//since no small selections, we save the big as the all
									roiManager("Deselect");
									RoiManager.translate(i*sizex, j*sizey);
									roiManager("Save", dir+image+"_" + i + "_"+ j +  "_ch" +channelld+"_all0.zip");
									selectWindow("ROI Manager");
									run("Close");
									
									
								} else {
									print("no adipocytes found, image must be wrong");
									error
								}
								
							}
						}
					}
				}
				
				//we have segmented all the crops, we now merge all segmentations in same image
				if (isOpen("ROI Manager")){
					selectWindow("ROI Manager");
					run("Close");
				}
				print("C");				
				selectWindow(image8);
				print("D");	
				
				if (ncropsx!=1 || ncropsy!=1 || dosmall!=false){
					for ( i=0; i<ncropsx; i++) {
						for ( j=0; j<ncropsy; j++) {
							if(File.exists(dir+image +"_" + i + "_"+ j +  "_ch" +channelld+ "_all0.zip")){
							roiManager("Open", dir+image +"_" + i + "_"+ j +  "_ch" +channelld+ "_all0.zip");
							}
						}
					}
					ntrois=roiManager("Count");
					if (ntrois>0){
					roiManager("Deselect");
					roiManager("Save", dir+image+  "_ch" +channelld+"_allcrops0.zip");
					}
				}
				if (isOpen("ROI Manager")){
					selectWindow("ROI Manager");
					run("Close");
				}
				run("Close");
				wait(500);
				
				if (measure==true){
					print("measuring");
					if(File.exists(dir+image+  "_ch" +channelld+"_allcrops0.zip")){
						roiManager("Open", dir+image+  "_ch" +channelld+"_allcrops0.zip");
						if(isOpen(image8p)){
							selectWindow(image8p);
						} else {
							print("when did we close image8p????");
							selectWindow(image);
							run("RGB Color");

							run("8-bit");
							image8p=getTitle();
						}

						roiManager("deselect");
						roiManager("Measure");
						selectWindow("Results");
						saveAs("results", dir+image+  "_ch" +channelld+"_RESULTSFINAL.txt");
						selectWindow("ROI Manager"); run("Close");
						selectWindow("Results"); run("Close");
						
					} else {
						print("We can't measure because, somehow, the final zip is missing");
					}
				}
				
				
				run("Close All");
			
			}
			
		}
	}
}
		
