macro "InDi3" {

//-----------------STARTING-----------------------
close("\\Others");
run("Close All");
// CLEANING RESULT TABLE
run("Clear Results");
// CLEANING LOG FILE
print("\\Clear");
// COLOR CODE
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);

//-----------------------------------------------

//----FILE LOCATION-----------------
open();
name_file = File.nameWithoutExtension;
full_name_file = File.name;
dir = getDirectory("Saving directory");
save_folder = dir + name_file + "_results/";
if ( File.exists(save_folder) < 1 ) File.makeDirectory( save_folder );
//----CLASSIFIER LOCATION-----------------
dir_classifier = "D://Nasim";
//dir_classifier = "G://PTEN/P24/Reelin-AKT/classifier" // copy-paste the classifier directory

//---INITIALIZATION-------------------
Stack.getDimensions(width, height, channels, slices, frames);
run("Channels Tool...");
Stack.setDisplayMode("composite");
run("Brightness/Contrast...");

Ch_nucleus = 1;
Ch_quantif1 = 2;
Ch_quantif2 = 3;
Ch_cell = 4;

temp_OG_hyperstack = getImageID(); run("Subtract Background...", "rolling=25");
Stack.setChannel(Ch_cell);
waitForUser("SUBSTACK DEFINITION","Remember Z_min and Z_max");
Z_min = getNumber("First Z", 1);
Z_max = getNumber("Last Z", slices);
run("Duplicate...", "duplicate slices="+Z_min+"-"+Z_max);
OG_hyperstack = getImageID();

selectImage(temp_OG_hyperstack); close();
selectImage(OG_hyperstack);




// --------------------------------------------
// --- 3D Cellular extraction
selectImage(OG_hyperstack);
run("Duplicate...", "duplicate channels=" +Ch_cell); // cell
stack_cell = getImageID();

//run("Gaussian Blur 3D...", "x=2 y=2 z=2");
run("Trainable Weka Segmentation 3D");
wait(1000);
//waitForUser("RESULT","LOAD AND CREATE");
call("trainableSegmentation.Weka_Segmentation.loadClassifier", dir_classifier+"/classifier.model");
wait(500);
call("trainableSegmentation.Weka_Segmentation.setFeature", "Gaussian_blur=true");
call("trainableSegmentation.Weka_Segmentation.setFeature", "Derivatives=true");
call("trainableSegmentation.Weka_Segmentation.setFeature", "Structure=true");
//call("trainableSegmentation.Weka_Segmentation.setFeature", "Difference_of_Gaussian=true");
call("trainableSegmentation.Weka_Segmentation.setFeature", "Minimum=true");
call("trainableSegmentation.Weka_Segmentation.setFeature", "Maximum=true");
call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
wait(1000);
call("trainableSegmentation.Weka_Segmentation.getResult");
close("Trainable*");
classified_cell = getImageID(); rename("Classified_cell");

selectImage(classified_cell);
run("Convert to Mask", "method=Huang background=Light calculate");
run("3D Fill Holes"); //run("Fill Holes", "stack");

setTool("Paintbrush Tool");
waitForUser("Manual Correction","Leave only cell bodies.\nPress Ok when you are done");


// -----------------
// --- 3D Section Contour extraction
setForegroundColor(255, 255, 255); setBackgroundColor(0, 0, 0);
selectImage(OG_hyperstack);
run("Duplicate...", "duplicate channels="+Ch_quantif2); // section contour
temp_projection_Contour = getImageID(); rename("Contour_section");
run("Subtract Background...", "rolling=50 disable stack");
run("Gaussian Blur 3D...", "x=2 y=2 z=2");

run("Convert to Mask", "method=Huang background=Dark calculate");
run("3D Fill Holes"); run("Fill Holes", "stack");
run("Invert", "stack"); run("Fill Holes", "stack"); run("Invert", "stack"); //cleaning the segmentation result

// --- Merging Resulting Extraction Files/Stacks
run("Concatenate...", "open image1=[Classified_cell] image2=[Contour_section] image3=[-- None --]");
run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
Stack.setDisplayMode("color"); Stack.setChannel(1); run("Grays"); Stack.setChannel(2); run("Grays");

rename("Cells_Contour"); saveAs("tiff", save_folder +"Cells_Contour");
Analysis_stack = getImageID();



//-----------------------------------------
//-----------------------------------------
//--------------- ROI EXTRACTION ----------

run("3D Manager"); Ext.Manager3D_Reset();
run("3D Manager Options", "volume surface compactness 3d_moments integrated_density mean_grey_value std_dev_grey_value minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest distance_between_centers=100 distance_max_contact=1.80 drawing=Contour use_0 use_1");

count_3dcells = 0;
count_3dcontour = 0;
//-----------------------------------------
//--- CELL 3D SEGMENTATION
selectImage(Analysis_stack);
Stack.setChannel(1); //Cells
Ext.Manager3D_Segment(128, 255); rename("Cell_3d_segment");
Ext.Manager3D_AddImage(); wait(2000); 
Seg_cells = getImageID();

		//--- USER QUALITY CONTROL FOR CELL 3D SEGMENTATION
selectImage(OG_hyperstack); Stack.setChannel(Ch_cell); run("Select None");
waitForUser("3D CELL SEGMENTATION","PRESS LIVE ON AND GO THROUGH ROIS"); run("Select None");
selectImage(Seg_cells); 
Stack.setChannel(1); run("Select None"); Ext.Manager3D_DeselectAll();
Ext.Manager3D_Count(count_3dcells);

if (count_3dcells < 1) {
	exit("You need at least one detected cell to proceed.\nMacro ends here.");
} else {
	selectImage(Seg_cells);
	for (i = 0; i < count_3dcells; i++) {
		Ext.Manager3D_Select(i);
		Ext.Manager3D_Rename("Cell_" +i);
		Ext.Manager3D_DeselectAll(); run("Select None");
	}
	Ext.Manager3D_Save(save_folder+"\Roi3D cell.zip"); wait(1000);
	Ext.Manager3D_DeselectAll(); Ext.Manager3D_Reset();
	
}
selectImage(OG_hyperstack); run("Select None");
selectImage(Seg_cells); run("Select None");
//-----------------------------------------
//--- CONTOUR 3D SEGMENTATION
selectImage(Analysis_stack); run("Select None"); 
run("3D Manager Options", "volume surface compactness 3d_moments integrated_density mean_grey_value std_dev_grey_value minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest distance_between_centers=100 distance_max_contact=1.80 drawing=Contour use_0 use_1");
Stack.setChannel(2); //Section contour
//run("Select None"); Ext.Manager3D_DeselectAll();
Ext.Manager3D_Segment(128, 255); rename("Contour_3d_segment");
Ext.Manager3D_AddImage(); wait(2000);
Seg_section = getImageID();

		//--- USER QUALITY CONTROL FOR CONTOUR 3D SEGMENTATION
selectImage(OG_hyperstack); Stack.setChannel(Ch_cell); run("Select None");
waitForUser("3D CONTOUR SEGMENTATION","PRESS LIVE ON AND GO THROUGH ROIS"); run("Select None");
selectImage(Seg_section);
Stack.setChannel(2); run("Select None"); Ext.Manager3D_DeselectAll();
Ext.Manager3D_Count(count_3dcontour);

if (count_3dcontour > 1) {
	waitForUser("3D CONTOUR SEGMENTATION","CHECK ONCE MORE");
}
Ext.Manager3D_Count(count_3dcontour);
if (count_3dcontour > 1) {
	exit("This macro only works with one section border!");
}
selectImage(Seg_section);
for (i = 0; i < count_3dcontour; i++){
	Ext.Manager3D_Select(i);
	Ext.Manager3D_Rename("Contour_" +i);
	Ext.Manager3D_DeselectAll(); run("Select None");
}
Ext.Manager3D_Save(save_folder+"\Roi3D contour.zip"); wait(2000);
Ext.Manager3D_DeselectAll(); Ext.Manager3D_Reset();

selectImage(OG_hyperstack); run("Select None");
selectImage(Seg_cells); run("Select None");

//-----------------------------------------
//--- CONCATENATION
selectImage(Analysis_stack); close();
run("Concatenate...", "open image1=[Cell_3d_segment] image2=[Contour_3d_segment] image3=[-- None --]");
run("Re-order Hyperstack ...", "channels=[Frames (t)] slices=[Slices (z)] frames=[Channels (c)]");
Stack.setDisplayMode("color"); Stack.setChannel(1); run("Grays"); Stack.setChannel(2); run("Grays");
rename("Cells_Contour_3d_segment"); saveAs("tiff", save_folder +"Cells_Contour_Seg");
Analysis_3d_stack = getImageID();

//-----------------------------------------
//-----------------------------------------
//--------------- QUANTIFICATION ----------
print("\\Clear"); run("Clear Results");
selectImage(Analysis_3d_stack); close();

Ext.Manager3D_Load(save_folder+"\Roi3D cell.zip"); // cell ROIs
wait(5000);
Ext.Manager3D_Count(count_3dcells); print(count_3dcells);
Ext.Manager3D_Load(save_folder+"\Roi3D contour.zip"); // contour ROI

selectImage(OG_hyperstack);

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// --- DISTANCE QUANTIFICATION
print("QUANTIFICATION");
count_all = 0; Ext.Manager3D_Count(count_all);// total number of ROIs
print(count_all);
dist_cellsToborder = newArray(count_3dcells); // array to collect distance
position_contour = count_3dcells; //position in 3D ROI MANAGER 
//waitForUser("BREAK","BREAK");
selectImage(OG_hyperstack); run("Select None");
//Ext.Manager3D_Distance(); // all distances are computed
//Ext.Manager3D_SaveResult("D",save_folder +"Distance.csv"); // all computed distances are saved
	
for (i=0; i <count_3dcells; i++){
		Ext.Manager3D_Dist2(i,count_3dcells,"c1b2",dist); //shortest distance from cell center to border
		dist_cellsToborder [i] = dist;
		Ext.Manager3D_DeselectAll();
	}
	
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// --- INTENSITY QUANTIFICATION -OVERALL
Ext.Manager3D_Reset();
selectImage(OG_hyperstack); run("Select None");

//---- CHANNEL 2------------------------------------------------------------
Ext.Manager3D_Load(save_folder+"\Roi3D cell.zip"); 
Stack.setChannel(2); //Reelin CHANNEL

Ext.Manager3D_Measure();
Ext.Manager3D_SaveResult("M",save_folder +"Ch2_ResultsShape.csv"); //only once for the shape measurements

Ext.Manager3D_Quantif();
Ext.Manager3D_SaveResult("Q",save_folder +"Ch2_ResultsIntensity.csv"); // Density, intensity ...

//---- CHANNEL 3------------------------------------------------------------
Ext.Manager3D_DeselectAll();
selectImage(OG_hyperstack); run("Select None");
Stack.setChannel(3);
Ext.Manager3D_Quantif();
Ext.Manager3D_SaveResult("Q",save_folder +"Ch3_ResultsIntensity.csv");
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// --- INTENSITY QUANTIFICATION -SUMMARY

ch2_intDen_3d_cell = newArray(count_3dcells);
ch3_intDen_3d_cell = newArray(count_3dcells);
ch2_mean_3d_cell = newArray(count_3dcells);
ch3_mean_3d_cell = newArray(count_3dcells);
ratio_mean = newArray(count_3dcells);
ratio_intDen = newArray(count_3dcells);
label = newArray(count_3dcells);
volume_3d_cell = newArray(count_3dcells);
surface_3d_cell = newArray(count_3dcells);
sphericity_3d_cell = newArray(count_3dcells);

selectImage(OG_hyperstack);

// --- Collecting results
for (i=0; i < count_3dcells; i++) {
	object = 0;
	object = object +i;
	// Extracting intDen and mean
	Ext.Manager3D_Select(object);
	Stack.setChannel(2);
	Ext.Manager3D_Quantif3D(object,"IntDen",intDen_cell_ch2);
	Stack.setChannel(2);
	Ext.Manager3D_Quantif3D(object,"Mean",mean_cell_ch2);
	ch2_intDen_3d_cell [object] = intDen_cell_ch2;
	ch2_mean_3d_cell [object] = mean_cell_ch2;

	Ext.Manager3D_GetName(object, name);
	label [object] = name;
	
	Ext.Manager3D_DeselectAll();
	Stack.setChannel(2); run("Select None");
	// Extracting intDen and mean
	Ext.Manager3D_Select(object);//
	Stack.setChannel(3);
	Ext.Manager3D_Quantif3D(object,"IntDen",intDen_cell_ch3);
	Stack.setChannel(3);
	Ext.Manager3D_Quantif3D(object,"Mean",mean_cell_ch3);
	ch3_intDen_3d_cell [object] = intDen_cell_ch3;
	ch3_mean_3d_cell [object] = mean_cell_ch3;
	Ext.Manager3D_DeselectAll();
	Stack.setChannel(3); run("Select None");
	// Extracting the volume, the surface, 
	Ext.Manager3D_Select(object);//
	Ext.Manager3D_Measure3D(object,"Vol",measure_volume);
	Ext.Manager3D_Measure3D(object,"Surf",measure_surface);
	Ext.Manager3D_Measure3D(object,"Spher",measure_sphericity);
	volume_3d_cell [object] = measure_volume;
	surface_3d_cell [object] = measure_surface;
	sphericity_3d_cell [object] = measure_sphericity;
//---------------------------------------------
// --- Building the summary table
//---------------------------------------------
	Ext.Manager3D_DeselectAll();
	//setResult("Name",object,"Cell_"+object);
	setResult("Name",object,label [object]);
	setResult("Ch2_IntDen3D",object,ch2_intDen_3d_cell [object]);
	setResult("Ch3_IntDen3D",object,ch3_intDen_3d_cell [object]);
	setResult("Ch2_Mean3D",object,ch2_mean_3d_cell [object]);
	setResult("Ch3_Mean3D",object,ch3_mean_3d_cell [object]);
	updateResults();
	
	ratio_intDen [object] = ch2_intDen_3d_cell [object] / ch3_intDen_3d_cell [object];
	ratio_mean [object] = ch2_mean_3d_cell [object] / ch3_mean_3d_cell [object];
	setResult("Ratio IntDen",object,ratio_intDen [object]);
	setResult("Ratio Mean",object,ratio_mean [object]);
	setResult("3d_Volume (um3)",object, volume_3d_cell [object]);
	setResult("3d_Surface (um2)", object, surface_3d_cell [object]);
	setResult("Sphericity", object, sphericity_3d_cell [object]);

	setResult("3D Dist to BORDER", object, dist_cellsToborder [object]);
	
	updateResults();
}
saveAs("results", save_folder +"Main Results.csv");

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// --- SURROUNDING INTENSITY QUANTIFICATION -REELIN CHANNEL
Ext.Manager3D_DeselectAll(); Ext.Manager3D_Reset();
run("Clear Results");
if (isOpen("ROI Manager")) { 
	roiManager("reset");
}

selectImage(OG_hyperstack);
getVoxelSize(width, height, depth, unit);

Ch_quantif1 = 2;
pixel_width = width; pixel_height = height; pixel_depth = depth;
voxel = pixel_width*pixel_width*pixel_depth;
inv_pixel_size = 1 / pixel_width;


run("Duplicate...", "duplicate channels=" +Ch_quantif1);
ref_stack = getImageID();
Stack.getDimensions(width, height, channels, slices, frames);
ref_width = width; ref_height = height;
ref_slices = slices;

selectImage(OG_hyperstack); close();
run("3D Manager");

Ext.Manager3D_Load(save_folder+"\Roi3D cell.zip"); wait(1000);
Ext.Manager3D_Count(count_3dcells);

factor = 5; // modulo *1um 
enlarging_factor = 0;
n_2Droi = 0; // counting the number of 2D ROIs
count = -1;

inital_V = newArray(count_3dcells);
effective_int = newArray(factor*count_3dcells);
effective_vol = newArray(factor*count_3dcells);

// ---------------------------------------------------------
// ---- Converting 3D to 2D ROI ----------------------------
for (i = 0; i < count_3dcells; i++) {
	newImage("Ref_3Dto2DROI_"+i, "8-bit black", ref_width, ref_height, ref_slices);
	Ext.Manager3D_Select(i);
	Ext.Manager3D_FillStack(255, 255, 255);
	Ext.Manager3D_DeselectAll(); run("Select None");
	run("Convert to Mask", "method=Default background=Dark calculate");
	run("Analyze Particles...", "clear add stack");	
	n_2Droi = roiManager("count");
	if (n_2Droi < 1) { exit("Not enough 2D ROI to continue"); } // WARNING
	roiManager("save", save_folder + "2DROI_Cell_#"+i+".zip");
	close("Ref_3Dto2DROI_"+i);

	temp_initial_V = 0;
	// ---- Collecting initial volume
	selectImage(ref_stack);
	for (k = 0; k < n_2Droi; k++) {
		roiManager("select",k);
		run("Clear", "slice");
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		temp_initial_V = temp_initial_V + nPixels;
	}
	roiManager("deselect");
	inital_V [i] = temp_initial_V*voxel; //in um	

	temp_intensity_V = 0;
	temp_volume = 0;
	// ---- Building donuts
	for (j = 1; j <= factor ; j++) {
		
		enlarging_factor = j * inv_pixel_size;
		count = count +1;
		
		for (k = 0; k < n_2Droi; k++) {
			roiManager("select",k);
			run("Enlarge...", "enlarge="+enlarging_factor);
			roiManager("update");
		}
		roiManager("save", save_folder + "2DROI_Cell_#"+i+"_Around_"+j+"um.zip");

		roiManager("Show None"); roiManager("deselect"); roiManager("Delete");
		roiManager("open", save_folder + "2DROI_Cell_#"+i+"_Around_"+j+"um.zip");
		// ---- Collecting DONUT intensity and volume
		for (k = 0; k < n_2Droi; k++) {
			roiManager("select",k);
			getRawStatistics(nPixels, mean, min, max, std, histogram);
			temp_volume = temp_volume + nPixels*voxel;
			temp_intensity_V = temp_intensity_V + nPixels*mean;
		}
		effective_vol [count] = temp_volume;
		effective_int [count] = temp_intensity_V;

	}
		
	roiManager("Show None"); roiManager("deselect"); roiManager("Delete");

}

// --------------- BUILDING RESULT TABLE ---------------------
// -----------------------------------------------------------

count = -1;

run("Clear Results");
for (n = 0; n < count_3dcells; n++) {
	for (m = 1; m <= factor ; m++) {
	count = count +1;
	setResult("Label", count, "Cell_#"+n+"_Around_"+m);
	setResult("Intensity", count, effective_int[count]);
	setResult("Volume",count, effective_vol[count]-inital_V [n]);
	setResult("Norm.Int", count, effective_int[count]/effective_vol[count]);	
	updateResults();
	}
	updateResults();
}
saveAs("results", save_folder +"REELIN secretion Results.csv");
//--------------------
selectImage(ref_stack); run("Select None");
saveAs("tiff", save_folder +"Reelin_Channel");
print("\\Clear");
}