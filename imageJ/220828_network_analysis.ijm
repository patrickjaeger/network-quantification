/* Patrick Jaeger; patrick.jaeger@hest.ethz.ch
 * Laboratory for Orthopaedic Biomechanics, ETH Zurich, Switzerland
 * Fiji/ImageJ 2.1.0/1.53c; Java 1.8.0_172 [64-bit]
 * 2022-08
 */

#@ File (label = "Input files", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "Dataset name", value = "network1") dataset
#@ String (label = "Image suffix", value = ".tif") suffix
#@ boolean (label = "Enhance images") enhance

// MAIN -----------------------------------------
run("Bio-Formats Macro Extensions");
out_dirs = createOutDirs(output, dataset);
processFolder(input);
save_nucleus_count(out_dirs[0], dataset);
print("FINISHED");

// FUNCTIONS ------------------------------------

function createOutDirs(input, dataset) {
  // Create output directories if they do not already exist
  resOutDir = input + File.separator + dataset + "_results_csv";
  imgOutDir = input + File.separator + dataset + "_results_img";
  
  if (!File.exists(input + File.separator + dataset + "_results_csv")) {
    File.makeDirectory(resOutDir);
    File.makeDirectory(imgOutDir);
  }
  
  dirs = newArray(2);
  dirs[0] = resOutDir;
  dirs[1] = imgOutDir;
  
  return dirs;
}

function processFolder(input) {
 list = getFileList(input);
  list = Array.sort(list);
  for (i = 0; i < list.length; i++) {
    if(File.isDirectory(input + File.separator + list[i]))
      processFolder(input + File.separator + list[i]);
    if(endsWith(list[i], suffix))
      processFile(input, list[i]);
  }
}

function processFile(input, file) { 
  BFopen(input, file);
  source_img = getTitle();
  max_proj_plus(source_img);
  if (enhance) {
    create_enhanced_composite(source_img);
  } else {
    create_raw_composite(source_img);
  }
  count_nuclei(source_img);
  quant_networks(source_img, out_dirs[0]);
  draw_network(source_img, out_dirs[1], dataset);
  close(source_img);
}

function BFopen(input, file) { 
  // Open input using the bioformats importer
  run("Bio-Formats Importer", 
  "open=[" + input + File.separator + file + 
  "] autoscale color_mode=Default rois_import=[ROI manager]" +
  " view=Hyperstack stack_order=XYCZT");
}

function max_proj_plus(target_img) { 
// Z-projection plus isolation of relevant channels
  selectWindow(target_img);
  run("Arrange Channels...", "new=12");  // C1: nuclei; C2: actin; C3: trash
  run("Z Project...", "projection=[Max Intensity]");
}

function selectPattern(pattern) {
  // Select window that contains the pattern anywhere 
  // in the title; no * required
  images = getList("image.titles");

  for (i = 0; i < images.length; i++){
    if (matches(images[i], ".*" + pattern + ".*")) selectWindow(images[i]);
  }
}

function create_enhanced_composite(target_img) { 
// Do some image enhancement and create a composite of the source image
  max_proj_plus(target_img);
  run("Split Channels");

  //// NUCLEI
  selectPattern("C1-");
  nucleus_img = getTitle();

  run("Median...", "radius=2");
  //run("Minimum...", "radius=2");
  //run("Maximum...", "radius=4");

  run("Bandpass Filter...", "filter_large=40 filter_small=15 suppress=None tolerance=5 autoscale saturate");
  run("Subtract Background...", "rolling=25");
  run("Enhance Contrast", "saturated=0.35");

  //// ACTIN
  selectPattern("C2-");
  actin_img = getTitle();
  
  run("Bandpass Filter...", "filter_large=200 filter_small=3 suppress=None tolerance=5 autoscale saturate");
  run("Subtract Background...", "rolling=25");
  run("Enhance Contrast", "saturated=0.35");
  run("Unsharp Mask...", "radius=2 mask=0.60");

  // Merge
  run("Merge Channels...", "c2=[" + actin_img + "] c3=[" + nucleus_img + "] create");
  run("Scale Bar...", "width=100 height=25 font=38 color=White background=None location=[Lower Right] bold overlay");
  run("RGB Color");
  close("Composite");
  close("MAX_" + target_img);
}
//create_enhanced_composite(getTitle());

function create_raw_composite(target_img) { 
// function description
  max_proj_plus(target_img);
  run("Split Channels");

  selectPattern("C1-");
  nucleus_img = getTitle();
  run("Enhance Contrast", "saturated=0.35");

  selectPattern("C2-");
  actin_img = getTitle();
  run("Enhance Contrast", "saturated=0.35");

  run("Merge Channels...", "c2=[" + actin_img + "] c3=[" + nucleus_img + "] create");
  run("Scale Bar...", "width=100 height=25 font=38 color=White background=None location=[Lower Right] bold overlay");
  run("RGB Color");
  close("Composite");
  close("MAX_" + target_img);
}
//create_raw_composite(getTitle());

function count_nuclei(target_img) { 
// Count nuclei
  max_proj_plus(target_img);
  run("Split Channels");
  close("C2-*");
  selectPattern("C1-");
  nucleus_img = getTitle();

  run("Median...", "radius=2");
  //run("Maximum...", "radius=4");
  run("Bandpass Filter...", "filter_large=40 filter_small=15 suppress=None tolerance=5 autoscale saturate");
  //run("Subtract Background...", "rolling=20 disable");
  
  setAutoThreshold("Otsu dark");
  run("Convert to Mask");
  run("Open");
  run("Watershed");

  run("Analyze Particles...", "size=30-Infinity circularity=0.50-1.00 " +
  "show=Overlay summarize add");

  // Draw outlines on composite image
  selectWindow("Composite (RGB)");
  roiManager("Show All without labels");
  setForegroundColor(255, 255, 0);  // yellow
  run("Line Width...", "line=1");
  roiManager("Draw");

  close(nucleus_img);
  roiManager("reset");
}

function quant_networks(target_img, output) { 
// Segment the actin network, then skeletonize and measure
  // Isolate channel
  max_proj_plus(target_img);
  run("Split Channels");
  close("C1-*");
  selectPattern("C2-");
  actin_img = getTitle();

  // Segmentation
  
  //// This was fine for study 1, but struggled with "particle noise" in study 2
  //run("Bandpass Filter...", "filter_large=50 filter_small=8 suppress=None tolerance=5 autoscale saturate");
  //setAutoThreshold("Huang dark");
  //setOption("BlackBackground", true);
  //run("Convert to Mask");

  //// This worked well in study 2, but is too sensitive for some images in study 1
  run("Bandpass Filter...", "filter_large=50 filter_small=7 suppress=None tolerance=5 autoscale saturate");
  run("Bandpass Filter...", "filter_large=50 filter_small=3 suppress=None tolerance=5 autoscale saturate");
  run("Subtract Background...", "rolling=25");
  run("Maximum...", "radius=5");
  setAutoThreshold("Huang dark");
  setOption("BlackBackground", true);
  run("Convert to Mask");

  // Remove some noise and skeletonize
  run("Close-");
  run("Skeletonize");

  // Remove length<10 pieces
  //run("Analyze Particles...", "size=10-Infinity pixel show=Masks");
  //run("Convert to Mask");
  //for (i = 0; i < 5; i++) run("Dilate");

  run("Analyze Skeleton (2D/3D)", "prune=none prune_0 show");
  /* Analyze Skeleton
   *  Output: 
   *  - Branch information
   *  - Results: - each row represents one connected network
   *             - sum junction voxels and slab voxels to get
   *               the total length of the network
   *  - Tagged skeleton: - shows all networks
   *                     - orange: branch pixels
   *                     - magenta: junctions
   *                     - purple: endpoints
   *                     - flatten, 8bit, make binary, 
   *                       draw ontop of source image to show network
   *  - Labeled skeletons: - assign each connnected network a color
   *                       - not really useful ATM; currently disabled
   */

  // Save branch and network information
  selectWindow("Branch information");
  branch_file_name = substring(target_img, 0, lengthOf(target_img)-4) + "_branch-info.csv";
  saveAs("Results", output + File.separator + branch_file_name);
  close(branch_file_name);

  selectWindow("Results");
  network_file_name = substring(target_img, 0, lengthOf(target_img)-4) + "_network-info.csv";
  saveAs("Results", output + File.separator + network_file_name);
  close("Results");

  close("C2-*");
}

function draw_network(target_img, output, out_name) { 
// Draw network on top of composite image and save
  selectWindow("Tagged skeleton");

  run("Flatten");
  run("8-bit");
  run("Make Binary");
  close("Tagged skeleton");
  selectWindow("Tagged skeleton-1");
  run("Create Selection");
  selectWindow("Composite (RGB)");
  setForegroundColor(255, 0, 255);  // magenta
  run("Line Width...", "line=1");
  run("Restore Selection");
  run("Draw", "slice");
  close("Tagged skeleton-1");
  
  selectWindow("Composite (RGB)");
  img_name = substring(target_img, 0, lengthOf(target_img)-4) + "_QQ.jpeg";
  saveAs("Jpeg", output + File.separator + img_name);
  close("Composite (RGB)");
}

function save_nucleus_count(output, dataset_name) { 
  selectWindow("Summary");
  table_name = dataset_name + "_nucleus_count.csv";
  saveAs("Results", output + File.separator + table_name);
  close(table_name);
}













