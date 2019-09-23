// provides
goog.provide('X.dicomUtils');

// requires
goog.require('X.base');

X.dicomUtils = {};

X.dicomUtils.BitsAllocated = Object.freeze({ "B8": 8, "B16": 16, "B24": 24, "B32": 32 });
X.dicomUtils.VOIFunction = Object.freeze({ "LINEAR": 0, "SIGMOID": 1 });
X.dicomUtils.PresentationLutShape = Object.freeze({ "IDENTITY": 0, "INVERSE": 1, "UNKNOWN": 2 });
X.dicomUtils.PhotometricInterpretation = Object.freeze({
    "MONOCHROME1": 0,
    "MONOCHROME2": 1,
    "RGB": 2,
    "PALETTECOLOR": 3,
    "UNKNOWN": 8
});
X.dicomUtils.PixelRepresentation = Object.freeze({ "UNSIGNED": 0, "SIGNED": 1 });

X.dicomUtils.ImageType = Object.freeze({ "JPEG": "jpeg", "DCM": "dcm", "MOBILE": "mobile", "CUSTOM": "custom" });

X.dicomUtils.applyPixelPaddingToHULookupTable = function (img, HULookupTable) {
    if (img.photometricInterpretation == X.dicomUtils.PhotometricInterpretation.MONOCHROME1 ||
        img.photometricInterpretation == X.dicomUtils.PhotometricInterpretation.MONOCHROME2) {

        if (img.pixelPaddingValue == null) {
            return;
        }

        var paddingValueMin = (img.pixelPaddingRange == null) ? img.pixelPaddingValue : Math.min(img.pixelPaddingValue, img.pixelPaddingRange);
        var paddingValueMax = (img.pixelPaddingRange == null) ? img.pixelPaddingValue : Math.max(img.pixelPaddingValue, img.pixelPaddingRange);

        var numPaddingValues = paddingValueMax - paddingValueMin + 1;
        var paddingValuesStartIndex = paddingValueMin;

        if (paddingValuesStartIndex >= HULookupTable.length) {
            return;
        }
        var isInverse = img.presentationLutShape == X.dicomUtils.PresentationLutShape.INVERSE ? true : img.photometricInterpretation == X.dicomUtils.PhotometricInterpretation.MONOCHROME1 ? true : false;
        var fill;
        if (isInverse) {
            fill = HULookupTable[HULookupTable.length - 1];
        } else {
            fill = HULookupTable[0];
        }

        for (var i = paddingValuesStartIndex; i < HULookupTable.length && i < paddingValuesStartIndex + numPaddingValues; i++) {
            HULookupTable[i] = fill;
        }
    }
};

X.dicomUtils.calculateHULookupTable = function (img) {

    var rs = img.rescale.getSlope();
    var ri = img.rescale.getIntercept();

    var bits = img.pixelPaddingValue != null ? img.bits : (rs == 1.0 && ri == 0.0 && img.pixelPaddingValue == null) ? img.bits : img.bitsStored;

    var size = Math.pow(2, bits);
    var HULookupTable = new Array(size);
    for (var i = 0; i < size; i++) {
        HULookupTable[i] = i * rs + ri;
    }

    X.dicomUtils.applyPixelPaddingToHULookupTable(img, HULookupTable);
    return HULookupTable;
};

X.dicomUtils.calculateLookupTable = function (img) {

    var HULookupTable = X.dicomUtils.calculateHULookupTable(img);

    var rs = img.rescale.getSlope();
    var ri = img.rescale.getIntercept();

    var size = Math.pow(2, img.pixelPaddingValue != null ? img.bits : (rs == 1.0 && ri == 0.0 && img.pixelPaddingValue == null) ? img.bits : img.bitsStored);

    var wc = img.windowLevel.getCenter();
    var ww = img.windowLevel.getWidth();

    // Create LUT it it does not exist
    var lookupTable = new Uint8Array(size);

    if (img.voiFunction == X.dicomUtils.VOIFunction.LINEAR) {

        var unsigned = img.pixelRepresentation == X.dicomUtils.PixelRepresentation.UNSIGNED;
        // Get offset and window limits
        var offset = (unsigned == true) ? 0 : size / 2,
            xMin = offset + wc - 0.5 - (ww - 1) / 2,
            xMax = offset + wc - 0.5 + (ww - 1) / 2,
            yMax = 255,
            yMin = 0;

        var minPixelValue = 0;
        var maxPixelValue = size;

        for (var inputValue = minPixelValue; inputValue < maxPixelValue; inputValue++) {
            if (HULookupTable[inputValue] <= xMin) {
                lookupTable[inputValue] = yMin;
            }
            else if (HULookupTable[inputValue] > xMax) {
                lookupTable[inputValue] = yMax;
            }
            else {
                var y = (((HULookupTable[inputValue] - offset) - (wc - 0.5)) / (ww - 1) + 0.5) * (yMax - yMin) + yMin;
                lookupTable[inputValue] = parseInt(y, 10);
            }
        }

    } else if (img.voiFunction == X.dicomUtils.VOIFunction.SIGMOID) {

        var outRangeSize = (1 << X.dicomUtils.BitsAllocated.B8) - 1;
        var maxOutValue = img.pixelRepresentation == X.dicomUtils.PixelRepresentation.SIGNED ? (1 << (X.dicomUtils.BitsAllocated.B8 - 1)) - 1 : outRangeSize;
        var minOutValue = img.pixelRepresentation == X.dicomUtils.PixelRepresentation.SIGNED ? -(maxOutValue + 1) : 0;
        var minInValue = 0;

        var nFactor = -20.0; // factor defined by default in Dicom standard ( -20*2/10 = -4 )
        var outRange = maxOutValue - minOutValue;
        var normalize = false;
        var minValue = 0, maxValue = 0, outRescaleRatio = 1;
        if (normalize) {
            var lowLevel = wc - ww / 2.0;
            var highLevel = wc + ww / 2.0;
            minValue = minOutValue + outRange / (1.0 + Math.exp((2.0 * nFactor / 10.0) * (lowLevel - wc) / ww));
            maxValue = minOutValue + outRange / (1.0 + Math.exp((2.0 * nFactor / 10.0) * (highLevel - wc) / ww));
            outRescaleRatio = (maxOutValue - minOutValue) / Math.abs(maxValue - minValue);
        }

        for (var i = 0; i < size; i++) {
            var value = outRange / (1.0 + Math.exp((2.0 * nFactor / 10.0) * (HULookupTable[i] + minInValue - wc) / ww));

            if (normalize) {
                value = (value - minValue) * outRescaleRatio;
            }

            value = parseInt(Math.round(value + minOutValue));
            value = parseInt((value > maxOutValue) ? maxOutValue : ((value < minOutValue) ? minOutValue : value));
            //value = parseInt(inverse ? (maxOutValue + minOutValue - value) : value);

            lookupTable[i] = value;
        }

    } else {
        throw ("Image VOI LUT not supported");
    }

    return lookupTable;
};

X.dicomUtils.initColorTableFromLUT = function (img, colorTable, lut) {

    var rs = img.rescale.getSlope();
    var ri = img.rescale.getIntercept();

    var offset;
    if (img.pixelRepresentation == X.dicomUtils.PixelRepresentation.UNSIGNED) {
        offset = 0;
    } else {
        offset = Math.pow(2, img.pixelPaddingValue != null ? img.bits : (rs == 1.0 && ri == 0.0 && img.pixelPaddingValue == null) ? img.bits : img.bitsStored) / 2;
    }

    // console.time('CALCULATE PIXELS');
    switch (img.photometricInterpretation) {
        case X.dicomUtils.PhotometricInterpretation.MONOCHROME1:
            for (var i = 0; i < lut.length; i++) {
                r = g = b = (255 - lut[i]);
                colorTable.add(i, 'lut', r, g, b, 255);
            }
            break;

        case X.dicomUtils.PhotometricInterpretation.MONOCHROME2:
            for (var i = 0; i < lut.length; i++) {
                r = g = b = lut[i];
                colorTable.add(i - offset, 'lut', r, g, b, 255);
            }
            break;
        default:
            throw new Error('X.dicomUtils.initColorTableFromLUT: photometricNotSupported');
    }
};

// export symbols (required for advanced compilation)
goog.exportSymbol('X.dicomUtils', X.dicomUtils);