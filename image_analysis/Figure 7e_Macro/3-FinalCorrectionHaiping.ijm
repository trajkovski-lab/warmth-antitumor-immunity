//This macro applies additional segmentation corrections based on area and circularity

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//minimum accepted values for mean intensity and circularity, and max accepted area
minmean=100; //was 160
maxarea=400;
mincirc=0.4;
//parameters for abhorrent shapes: segmentations that are bigger and less circular than these will be removed
maxareac=300;
mincircc=0.6;
//////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		


setOption("ExpandableArrays", true);
run("Set Measurements...", "area mean standard min centroid perimeter bounding fit shape feret's integrated skewness area_fraction display redirect=None decimal=5");

dir = getDirectory("Select a directory containing one or several lif files.");
files = getFileList(dir);

resolution=3;
seriesName = "all";


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

			finalname=file+ " - " + filenaming + " " + iseriesname + "_FinalResults.txt";
			print(finalname);
			
			if(!File.exists(finalname)) {
				run("Bio-Formats Importer", "open=[" + file + "] autoscale color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_list=" + seriesToOpen[s]);
				imageName=getTitle();
				title=imageName;
				run("RGB Color");
				run("8-bit");
				title8=getTitle();
				
				getDimensions(widthc, heightc, channelsc, slicesc, framesc);
				
				newImage("Untitled", "8-bit white", widthc, heightc, 1); 
				
				//we correct here the unwanted areas
				if (File.exists(dir + title + "_toeliminate.zip")){
					roiManager("Open", dir + title + "_toeliminate.zip");
					selectWindow(title8);
					setForegroundColor(0, 0, 0);
					setBackgroundColor(255, 255, 255);
					roiManager("deselect");
					roiManager("Fill");
					selectWindow("ROI Manager");
					run("Close");
				}
				



				roiManager("Open", dir+ title+ "__allcrops1.zip");
				
				print("Fiding and eliminating dark or huge selections");

				selectWindow(title8);
				roiManager("Deselect");
				roiManager("Measure");
				errors=newArray(0);
				totalrois=nResults();
				for (i=0; i<totalrois; i++) {

					mean=getResult("Mean",i);
					area=getResult("Area", i);
					circ = getResult("Circ.", i);

					sol=getResult("Solidity", i);
					widthi=getResult("Width", i);
					heighti=getResult("Height", i);
					

					areasq=widthi*heighti;
					
					if( sol<0.75 || (sol>0.99 && area > 0.99*areasq)  ||  (mean <minmean ) || (area> maxarea) || (area > maxareac && circ < mincircc) || (circ < mincirc)   ) { 
						errors= Array.concat(errors, i);
						
						if( (sol>0.99 && area > 0.99*areasq) ) { 

							roiManager("Select", i);
							roiManager("Set Color", "red");
							roiManager("Set Line Width", 2);

						}
					} else { 
						if ( area < 2 && mean > 195) { 
							selectWindow(title8);
							roiManager("Select", i);
							run("Fit Rectangle");
							run("Measure");
							newnres=nResults();
							newnres=newnres-1;
							areaR=getResult("Area", newnres);

							if (area/areaR > 0.9) {
								errors= Array.concat(errors, i);
								roiManager("Select", i);
								roiManager("Set Color", "cyan");
								roiManager("Set Line Width", 2);
							}
						}
					} 
				}

				
				
				nerrors=errors.length;
				if(nerrors>0){
					roiManager("deselect");
					roiManager("select", errors);
					roiManager("delete");
				}
				
				roiManager("deselect");
				selectWindow("Results");
				run("Close");
				
				roiManager("Show None");
				roiManager("Show All");
				setForegroundColor(0, 0, 0);
				setBackgroundColor(255, 255, 255);
				
				selectWindow("Untitled");
				roiManager("Fill");
				
				
	
				
				
				run("Distance Map");
				
				setThreshold(0, 3, "raw");
				run("Analyze Particles...", "size=250000-Infinity show=Masks clear include");
				roiManager("Show All");	
				
				selectWindow("Mask of Untitled");
				roiManager("Deselect");
				roiManager("Measure");
				
				
				print("Finding isolated selections");
				errors=newArray(0);
				for (i=0; i<nResults; i++) {
	
					mean=getResult("Mean",i);
					area=getResult("Area", i);
					if( (mean==0 ) || (area==0)   ) { 
						errors= Array.concat(errors, i);
					}
				}
				
				nerrors=errors.length;
				if(nerrors>0){
					roiManager("deselect");
					roiManager("select", errors);
					roiManager("delete");
				}
				selectWindow("Results");
				run("Close");
				roiManager("Deselect");
				
				
				nfrois=roiManager("Count");
				if(nfrois >0) {
					roiManager("save", dir+ title+"_all3.zip");
					selectWindow(title);
					roiManager("Measure");
					selectWindow("Results");
					saveAs("Results", dir+title+"_FinalResults.txt");
				}else {
					print("no adipocytes remaining in this image");
				}
				
			if (isOpen("ROI Manager")){
				selectWindow("ROI Manager");
				run("Close");
			}
			
			if (isOpen("Results")){
				selectWindow("Results");
				run("Close");
			}
			
			run("Close All");						
			}	
		}
	}
}
		


