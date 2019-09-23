/*
 * 
 *                  xxxxxxx      xxxxxxx
 *                   x:::::x    x:::::x 
 *                    x:::::x  x:::::x  
 *                     x:::::xx:::::x   
 *                      x::::::::::x    
 *                       x::::::::x     
 *                       x::::::::x     
 *                      x::::::::::x    
 *                     x:::::xx:::::x   
 *                    x:::::x  x:::::x  
 *                   x:::::x    x:::::x 
 *              THE xxxxxxx      xxxxxxx TOOLKIT
 *                    
 *                  http://www.goXTK.com
 *                   
 * Copyright (c) 2012 The X Toolkit Developers <dev@goXTK.com>
 *                   
 *    The X Toolkit (XTK) is licensed under the MIT License:
 *      http://www.opensource.org/licenses/mit-license.php
 * 
 *      'Free software' is a matter of liberty, not price.
 *      'Free' as in 'free speech', not as in 'free beer'.
 *                                         - Richard M. Stallman
 * 
 * 
 * CREDITS: Thank you to Thomas J. Re for his initial implementation.
 *
 */
// provides
goog.provide('X.parserIMAGEJS');
// requires
goog.require('X.event');
goog.require('X.object');
goog.require('X.parser');
goog.require('X.triplets');
goog.require('X.dicomUtils');
goog.require('goog.math.Vec3');
/**
 * Create a parser for DICOM files.
 * 
 * @constructor
 * @extends X.parser
 */
X.parserIMAGEJS = function () {
  //
  // call the standard constructor of X.parser
  goog.base(this);
  //
  // class attributes
  /**
   * @inheritDoc
   * @const
   */
  this._classname = 'parserIMAGEJS';
};
// inherit from X.parser
goog.inherits(X.parserIMAGEJS, X.parser);
/**
 * @inheritDoc
 */
X.parserIMAGEJS.prototype.parse = function (container, object, data, flag) {
  // X.TIMER(this._classname + '.parse');
  // needed, for renderer2d and 3d legacy...

  object.MRI = {};
  object.MRI.loaded_files = 0;

  // parse the byte stream
  this.parseStream(data, object);

  // return;
  // X.TIMERSTOP(this._classname + '.parse');
  // check if all slices were completed loaded
  if (!goog.isDefAndNotNull(object._file.length) || object.slices.length == object._file.length) {

    // needed, for renderer2d and 3d legacy...
    object.MRI.loaded_files = object._file.length;

    // sort slices per series
    var series = {};
    var imageSeriesPushed = {};
    for (var i = 0; i < object.slices.length; i++) {

      // series undefined yet
      if (!series.hasOwnProperty(object.slices[i]['series_instance_uid'])) {

        series[object.slices[i]['series_instance_uid']] = new Array();
        imageSeriesPushed[object.slices[i]['series_instance_uid']] = {};

      }

      // push image if it has not been pushed yet
      if (!imageSeriesPushed[object.slices[i]['series_instance_uid']].hasOwnProperty(object.slices[i]['sop_instance_uid'])) {

        imageSeriesPushed[object.slices[i]['series_instance_uid']][object.slices[i]['sop_instance_uid']] = true;
        series[object.slices[i]['series_instance_uid']].push(object.slices[i]);

      }

    }

    imageSeriesPushed = null;

    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    ////////////////////////////////////////////////////////////////////////

    // GLOBAL PARAMETERS
    // pointer to first image
    var seriesInstanceUID = Object.keys(series)[0];
    var first_image = series[seriesInstanceUID];
    // number of unique slices available
    var first_image_stacks = first_image.length;
    // container for volume specific information
    var volumeAttributes = {};

    ////////////////////////////////////////////////////////////////////////
    //
    // ORDER SLICES
    //
    ////////////////////////////////////////////////////////////////////////

    //
    // we can order slices based on
    //
    // image_position_patient:
    // -> each slice show have a different 'image_position_patient'
    // -> The Image Position (0020,0032) specifies the x, y, and z coordinates of
    // -> -> the upper left hand corner of the image; it is the center of the first
    // -> -> voxel transmitted. Image Orientation (0020,0037) specifies the direction
    // -> -> cosines of the first row and the first column with respect to the patient.
    // -> -> These Attributes shall be provide as a pair. Row value for the x, y, and
    // -> -> z axes respectively followed by the Column value for the x, y, and z axes
    // -> -> respectively.
    //
    // in some cases, such as diffusion, 'image_position_patient' is the same for all
    // slices. We should then use the instance_number to order the slices.
    // 
    // instance_number:
    // -> each slice show have a different 'instance_number'
    // -> A number that identifies this raw data. 
    // -> -> The value shall be unique within a series

    var _ordering = 'image_position_patient';


    if (first_image_stacks == 1) {

      // ORDERING BASED ON IMAGE POSITION
      _ordering = 'image_position_patient';

      // set distance to 0
      series[seriesInstanceUID][0]['dist'] = 0;

    }
    else if (first_image[0]['image_position_patient'][0] != first_image[1]['image_position_patient'][0] ||
      first_image[0]['image_position_patient'][1] != first_image[1]['image_position_patient'][1] ||
      first_image[0]['image_position_patient'][2] != first_image[1]['image_position_patient'][2]) {

      // ORDERING BASED ON IMAGE POSITION
      _ordering = 'image_position_patient';

      // set distances
      var _x_cosine = new goog.math.Vec3(first_image[0]['image_orientation_patient'][0],
        first_image[0]['image_orientation_patient'][1], first_image[0]['image_orientation_patient'][2]);
      var _y_cosine = new goog.math.Vec3(first_image[0]['image_orientation_patient'][3],
        first_image[0]['image_orientation_patient'][4], first_image[0]['image_orientation_patient'][5]);
      var _z_cosine = goog.math.Vec3.cross(_x_cosine, _y_cosine);

      function computeDistance(flag, arrelem) {
        arrelem['dist'] = arrelem['image_position_patient'][0] * flag.x +
          arrelem['image_position_patient'][1] * flag.y +
          arrelem['image_position_patient'][2] * flag.z;
        return arrelem;
      }

      // compute dist in this series
      first_image.map(computeDistance.bind(null, _z_cosine));
      // order by dist
      first_image.sort(function (a, b) { return a["dist"] - b["dist"] });

    }
    else if (first_image[0]['instance_number'] != first_image[1]['instance_number']) {

      // ORDERING BASED ON instance number
      _ordering = 'instance_number';
      first_image.sort(function (a, b) { return a["instance_number"] - b["instance_number"] });

    }
    else {

      window.console.log("Could not resolve the ordering mode");

    }


    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    ////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////
    //
    // COMPUTE SPACING
    //
    ////////////////////////////////////////////////////////////////////////

    if (isNaN(first_image[0]['pixel_spacing'][0])) {

      first_image[0]['pixel_spacing'][0] = 1.;

    }

    if (isNaN(first_image[0]['pixel_spacing'][1])) {

      first_image[0]['pixel_spacing'][1] = 1.;

    }

    if (first_image_stacks > 1) {

      switch (_ordering) {
        case 'image_position_patient':
          // We work only on 2 first slices
          var _first_position = first_image[0]['image_position_patient'];
          var _second_image_position = first_image[1]['image_position_patient'];
          var _x = _second_image_position[0] - _first_position[0];
          var _y = _second_image_position[1] - _first_position[1];
          var _z = _second_image_position[2] - _first_position[2];
          first_image[0]['pixel_spacing'][2] = Math.sqrt(_x * _x + _y * _y + _z * _z);
          break;
        case 'instance_number':
          first_image[0]['pixel_spacing'][2] = 1.0;
          break;
        default:
          window.console.log("Unkown ordering mode - returning: " + _ordering);
          break;
      }

    }
    else {

      first_image[0]['pixel_spacing'][2] = 1.0;

    }

    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    // -> we estimated the spacing in all directions
    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    //
    // Estimate number of slices we are expecting
    //
    ////////////////////////////////////////////////////////////////////////

    // we execpt at least one image :)
    var first_image_expected_nb_slices = 1;
    switch (_ordering) {
      case 'image_position_patient':
        // get distance between 2 points
        var _first_position = first_image[0]['image_position_patient'];
        var _last_image_position = first_image[first_image_stacks - 1]['image_position_patient'];
        var _x = _last_image_position[0] - _first_position[0];
        var _y = _last_image_position[1] - _first_position[1];
        var _z = _last_image_position[2] - _first_position[2];
        var _distance_position = Math.sqrt(_x * _x + _y * _y + _z * _z);
        //normalize by z spacing
        first_image_expected_nb_slices += Math.round(_distance_position / first_image[0]['pixel_spacing'][2]);
        break;
      case 'instance_number':
        first_image_expected_nb_slices += Math.abs(first_image[first_image_stacks - 1]['instance_number'] - first_image[0]['instance_number']);
        break;
      default:
        window.console.log("Unkown ordering mode - returning: " + _ordering);
        break;
    }

    var first_slice_size = first_image[0]['columns'] * first_image[0]['rows'];
    var first_image_size = first_slice_size * (first_image_expected_nb_slices);

    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    // -> we estimated the spacing in all directions
    // -> we know how many slices we expect in the best case
    ////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////
    //
    // Prepare and fill data container
    //
    ////////////////////////////////////////////////////////////////////////

    var first_image_data = null;

    // create data container

    if (first_image[0].pixel_representation) {
      switch (first_image[0].bits_allocated) {
        case 8:
          first_image_data = new Int8Array(first_image_size);
          break;
        case 16:
          first_image_data = new Int16Array(first_image_size);
          break;
        case 32:
          first_image_data = new Int32Array(first_image_size);
        default:
          window.console.log("Unknown number of bits allocated - using default: 32 bits");
          break;
      }
    } else {
      switch (first_image[0].bits_allocated) {
        case 8:
          first_image_data = new Uint8Array(first_image_size);
          break;
        case 16:
          first_image_data = new Uint16Array(first_image_size);
          break;
        case 32:
          first_image_data = new Uint32Array(first_image_size);
        default:
          window.console.log("Unknown number of bits allocated - using default: 32 bits");
          break;
      }
    }

    object._spacing = first_image[0]['pixel_spacing'];

    // fill data container
    // by pushing slices where we expect them
    // 
    // for instance, we have 3 non-consecutive slices
    // we are expecting 4 slices, the 3rd one is missing
    //
    // BEFORE:
    //
    // 0000000
    // 0000000
    // 0000000
    // 0000000
    //
    // AFTER:
    //
    // 1234123 -> first slice
    // 1234211 -> second slice
    // 0000000
    // 1232414 -> third slice

    for (var _i = 0; _i < first_image_stacks; _i++) {
      // get data
      var _data = first_image[_i].data;
      var _distance_position = 0;

      switch (_ordering) {
        case 'image_position_patient':
          var _x = first_image[_i]['image_position_patient'][0] - first_image[0]['image_position_patient'][0];
          var _y = first_image[_i]['image_position_patient'][1] - first_image[0]['image_position_patient'][1];
          var _z = first_image[_i]['image_position_patient'][2] - first_image[0]['image_position_patient'][2];
          _distance_position = Math.round(Math.sqrt(_x * _x + _y * _y + _z * _z) / first_image[0]['pixel_spacing'][2]);
          break;
        case 'instance_number':
          _distance_position = first_image[_i]['instance_number'] - first_image[0]['instance_number'];
          break;
        default:
          window.console.log("Unkown ordering mode - returning: " + _ordering);
          break;
      }

      first_image_data.set(_data, _distance_position * first_slice_size);
    }

    volumeAttributes.data = first_image_data;
    object._data = first_image_data;

    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    // -> we estimated the spacing in all directions
    // -> we know how many slices we expect in the best case
    // -> data container contains ordered data!
    ////////////////////////////////////////////////////////////////////////

    // IJK image dimensions
    // NOTE:
    // colums is index 0
    // rows is index 1
    object._dimensions = [first_image[0]['columns'], first_image[0]['rows'], first_image_expected_nb_slices];
    volumeAttributes.dimensions = object._dimensions;

    // get the min and max intensities
    var min_max = this.arrayMinMax(first_image_data);
    var min = min_max[0];
    var max = min_max[1];

    // attach the scalar range to the volume
    volumeAttributes.min = object._min = object._windowLow = min;
    volumeAttributes.max = object._max = object._windowHigh = max;
    // .. and set the default threshold
    // only if the threshold was not already set
    if (object._lowerThreshold == -Infinity) {

      object._lowerThreshold = min;

    }
    if (object._upperThreshold == Infinity) {

      object._upperThreshold = max;

    }

    // Slices are ordered so
    // volume origin is the first slice position
    var _origin = first_image[0]['image_position_patient'];

    //
    // Generate IJK To RAS matrix and other utilities.
    //

    var IJKToRAS = goog.vec.Mat4.createFloat32();

    ////////////////////////////////////////////////////////////////////////
    //
    // IMPORTANT NOTE:
    //
    // '-' added for LPS to RAS conversion
    // IJKToRAS is Identity if we have a time series
    //
    ////////////////////////////////////////////////////////////////////////

    if (object['reslicing'] == 'false' || object['reslicing'] == false) {
      goog.vec.Mat4.setRowValues(IJKToRAS,
        0,
        first_image[0]['pixel_spacing'][0],
        0,
        0,
        0);
      // - first_image[0]['pixel_spacing'][0]/2);
      goog.vec.Mat4.setRowValues(IJKToRAS,
        1,
        0,
        first_image[0]['pixel_spacing'][1],
        0,
        0);
      // - first_image[0]['pixel_spacing'][1]/2);
      goog.vec.Mat4.setRowValues(IJKToRAS,
        2,
        0,
        0,
        first_image[0]['pixel_spacing'][2],
        0);
      // + first_image[0]['pixel_spacing'][2]/2);
      goog.vec.Mat4.setRowValues(IJKToRAS,
        3, 0, 0, 0, 1);
    }
    else {
      switch (_ordering) {
        case 'image_position_patient':

          var _x_cosine = new goog.math.Vec3(first_image[0]['image_orientation_patient'][0],
            first_image[0]['image_orientation_patient'][1], first_image[0]['image_orientation_patient'][2]);
          var _y_cosine = new goog.math.Vec3(first_image[0]['image_orientation_patient'][3],
            first_image[0]['image_orientation_patient'][4], first_image[0]['image_orientation_patient'][5]);
          var _z_cosine = goog.math.Vec3.cross(_x_cosine, _y_cosine);

          goog.vec.Mat4.setRowValues(IJKToRAS,
            0,
            -first_image[0]['image_orientation_patient'][0] * first_image[0]['pixel_spacing'][0],
            -first_image[0]['image_orientation_patient'][3] * first_image[0]['pixel_spacing'][1],
            -_z_cosine.x * first_image[0]['pixel_spacing'][2],
            -_origin[0]);
          // - first_image[0]['pixel_spacing'][0]/2);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            1,
            -first_image[0]['image_orientation_patient'][1] * first_image[0]['pixel_spacing'][0],
            -first_image[0]['image_orientation_patient'][4] * first_image[0]['pixel_spacing'][1],
            -_z_cosine.y * first_image[0]['pixel_spacing'][2],
            -_origin[1]);
          // - first_image[0]['pixel_spacing'][1]/2);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            2,
            first_image[0]['image_orientation_patient'][2] * first_image[0]['pixel_spacing'][0],
            first_image[0]['image_orientation_patient'][5] * first_image[0]['pixel_spacing'][1],
            _z_cosine.z * first_image[0]['pixel_spacing'][2],
            _origin[2]);
          // + first_image[0]['pixel_spacing'][2]/2);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            3, 0, 0, 0, 1);
          break;
        case 'instance_number':
          goog.vec.Mat4.setRowValues(IJKToRAS,
            0, -1, 0, 0, -_origin[0]);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            1, -0, -1, -0, -_origin[1]);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            2, 0, 0, 1, _origin[2]);
          goog.vec.Mat4.setRowValues(IJKToRAS,
            3, 0, 0, 0, 1);
          break;
        default:
          window.console.log("Unkown ordering mode - returning: " + _ordering);
          break;
      }
    }

    volumeAttributes.IJKToRAS = IJKToRAS;
    volumeAttributes.RASToIJK = goog.vec.Mat4.createFloat32();
    goog.vec.Mat4.invert(volumeAttributes.IJKToRAS, volumeAttributes.RASToIJK);

    ////////////////////////////////////////////////////////////////////////
    // At this point:
    // -> slices are ordered by series
    // -> slices within a series are unique
    // -> we estimated the spacing in all directions
    // -> we know how many slices we expect in the best case
    // -> data container contains ordered data!
    // -> IJK To RAS (and invert) matrices
    ////////////////////////////////////////////////////////////////////////

    //
    // compute last required information for reslicing
    //

    // get RAS spacing
    //
    var tar = goog.vec.Vec4.createFloat32FromValues(0, 0, 0, 1);
    var res = goog.vec.Vec4.createFloat32();
    goog.vec.Mat4.multVec4(IJKToRAS, tar, res);

    var tar2 = goog.vec.Vec4.createFloat32FromValues(1, 1, 1, 1);
    var res2 = goog.vec.Vec4.createFloat32();
    goog.vec.Mat4.multVec4(IJKToRAS, tar2, res2);

    volumeAttributes.RASSpacing = [res2[0] - res[0], res2[1] - res[1], res2[2] - res[2]];

    // get RAS Boundung Box
    //
    var _rasBB = X.parser.computeRASBBox(IJKToRAS, [object._dimensions[0], object._dimensions[1], object._dimensions[2]]);
    // grab the RAS Dimensions
    volumeAttributes.RASDimensions = [_rasBB[1] - _rasBB[0] + 1, _rasBB[3] - _rasBB[2] + 1, _rasBB[5] - _rasBB[4] + 1];

    // get RAS Origin
    // (it is actually RAS min x, min y and min z)
    //
    volumeAttributes.RASOrigin = [_rasBB[0], _rasBB[2], _rasBB[4]];


    var centerImage = object.slices[Math.ceil(object.slices.length / 2)].ref_image;
    X.dicomUtils.initColorTableFromLUT(centerImage, object.colortable, X.dicomUtils.calculateLookupTable(centerImage));

    // create the volume object
    object.create_(volumeAttributes);

    for(var i = 0; i < object.slices.length; i++){
        delete object.slices[i]['ref_image'];
        delete object.slices[i]['data'];
    }
    delete object['slices'];
    volumeAttributes = null;
    first_image = null;
    first_image_data = null;
    series = null;

    // re-slice the data in SAGITTAL, CORONAL and AXIAL directions
    object._image = this.reslice(object);
    object._data = null;

  }

  // the object should be set up here, so let's fire a modified event
  var modifiedEvent = new X.event.ModifiedEvent();
  modifiedEvent._object = object;
  modifiedEvent._container = container;
  this.dispatchEvent(modifiedEvent);

};


/**
 * Parse the data stream according to the .nii/.nii.gz file format and return an
 * MRI structure which holds all parsed information.
 * 
 * @param {! PACScenter Image}
 *          data The PACScenter Image.
 * @param {!X.object}
 *          object The parent object.
 * @return {Object} The MRI structure which holds all parsed information.
 */
X.parserIMAGEJS.prototype.parseStream = function (data, object) {
  // attach the given data
  //object._data = data;

  var image = data.image;
  var pixelData = data.pixelData;

  if (typeof (object.slices) == "undefined" || object.slices == null) {
    object.slices = new Array();
  }

  // set slice default minimum required parameters
  var slice = {};
  slice['pixel_spacing'] = [.1, .1, Infinity];
  slice['image_orientation_patient'] = [1, 0, 0, 0, 1, 0];
  slice['image_position_patient'] = [0, 0, 0];
  // Transfer syntax UIDs
  // 1.2.840.10008.1.2: Implicit VR Little Endian
  // 1.2.840.10008.1.2.1: Explicit VR Little Endian
  // 1.2.840.10008.1.2.2: Explicit VT Big Endian
  slice['transfer_syntax_uid'] = "1.2.840.10008.1.2"; //TODO: check this value

  slice['rows'] = image.size.getHeight();
  slice['columns'] = image.size.getWidth();
  slice['bits_allocated'] = image.bits;
  slice['pixel_representation'] = image.pixelRepresentation;
  slice['photometric_interpretation'] = image.photometricInterpretation;
  slice['bits_stored'] = image.bitsStored;
  slice['number_of_images'] = image.series.images.length;
  slice['pixel_spacing'] = [image.spacing.getRowSpacing(), image.spacing.getColumnsSpacing(), Infinity];

  slice['series_instance_uid'] = image.series.id;
  slice['instance_number'] = image.instanceNumber;
  slice['image_position_patient'] = image.patientPosition.toArray();
  slice['image_orientation_patient'] = image.patientOrientation.toArray();

  slice['sop_instance_uid'] = image.id;

  if (image.pixelRepresentation == X.dicomUtils.PixelRepresentation.UNSIGNED) {
    slice['offset'] = 0;
  } else {
    var rs = image.rescale.getSlope();
    var ri = image.rescale.getIntercept();
    slice['offset'] = Math.pow(2, image.pixelPaddingValue != null ? image.bits : (rs == 1.0 && ri == 0.0 && image.pixelPaddingValue == null) ? image.bits : image.bitsStored) / 2;
  }

  // check for data type and parse accordingly
  var _data = null;
  ////Get pixels
  switch (slice.bits_allocated) {
    case X.dicomUtils.BitsAllocated.B8:
      if (image.pixelRepresentation === X.dicomUtils.PixelRepresentation.UNSIGNED) {
        _data = new Uint8Array(pixelData);
      } else {
        _data = new Int8Array(pixelData);
      }
      break;
    case X.dicomUtils.BitsAllocated.B16:
      if (image.pixelRepresentation === X.dicomUtils.PixelRepresentation.UNSIGNED) {
        _data = new Uint16Array(pixelData);
      } else {
        _data = new Int16Array(pixelData);
      }
      break;
    default:
      throw new Error('Allocated bits are currently not supported: ' + slice.bits_allocated);
  }

  slice['data'] = _data;
  slice['ref_image'] = data.image;

  object.slices.push(slice);

  return object;
};

// export symbols (required for advanced compilation)
goog.exportSymbol('X.parserIMAGEJS', X.parserIMAGEJS);
goog.exportSymbol('X.parserIMAGEJS.prototype.parse', X.parserIMAGEJS.prototype.parse);
