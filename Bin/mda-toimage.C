// ==========================================================================
// $Id: mda-toimage.C 253 2008-12-01 04:58:08Z heidrich $
// convert MDA stream to image using the ImageMagick library
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: heidrich (Wolfgang Heidrich)
// Email:   heidrich@cs.ubc.ca
// ==========================================================================

#include <Magick++.h>
#include <stdio.h>

#include "MDA/Base/ChannelList.hh"
#include "MDA/Array/MDAFileIO.hh"

using namespace MDA;
using namespace std;
using namespace Magick;


int
main(int argc, char *argv[])
{
	unsigned long i, j;

	CommandlineParser parser;
	
	// setup options
	ChannelList	channels;
	ChannelListOption channelOption(channels);
	parser.registerOption(&channelOption);
	
	// parse options
	int index = 1;
	if (!parser.parse(index, argc, argv) || index != argc - 1)
	{
		parser.usage(argv[0], "<options> <outfile>");
		exit(1);
	}

	// setup MDA reader and read input header
	MDAReader reader;
	reader.connect(cin);
	if (!reader.readHeader())
	{
		cerr << "Cannot read file header\n";
		exit(1);
	}
	
	DataType	type = reader.getType();
	CoordinateVector	dim = reader.getDim();
	unsigned int	numChannels = reader.getNumChannels();
	unsigned long numScanlines = reader.getNumScanlinesLeft();
	unsigned long scanlineSizeIn = reader.getScanlineSize();
	void *scanlineIn;

	unsigned int numOutChannels = channels.vec.size();
	if (numOutChannels == 0)
		// no channels specified? use the first few channels!
	for (numOutChannels = 0;
		numOutChannels< numChannels && numOutChannels< 4;
		numOutChannels++)
		channels.vec.push_back(numOutChannels);
	else
	{
		// check that we don't have too many channels
		if (numOutChannels> 4)
		{
			cerr << "Can only deal with 1-4 channels\n";
			exit(1);
		}
		// check that all channels in the channel list are vaild
		for (i = 0; i< numOutChannels; i++)
		if (channels.vec[i] >= numChannels)
		{
			cerr << "Input MDA only has " << numChannels << " channels\n";
			exit(1);
		}
	}

	// setup ImageMagick image
	ostringstream geometry;
	geometry << dim.vec[0] << 'x' << numScanlines;
	Image outImage(geometry.str().c_str(), "black");

	switch (numOutChannels)
	{
		case 1:
		outImage.type(GrayscaleType);
		break;
	case 2:
#if defined (_WIN32) || defined (_WIN64)
		outImage.type(GrayscaleAlphaType);
#else
		outImage.type(GrayscaleMatteType);
#endif
		break;
	case 3:
		outImage.type(TrueColorType);
		break;
	case 4:
#if defined (_WIN32) || defined (_WIN64)
		outImage.type(TrueColorAlphaType);
#else
		outImage.type(TrueColorMatteType);
#endif
		break;
	}
	outImage.modifyImage();
	Pixels view(outImage);


	unsigned short *scanlineOut = new unsigned short[dim.vec[0] * numChannels];

	// copy image data
	for (i = 0; i< numScanlines; i++)
	{
		// read scanline and convert to unsigned short
		scanlineIn = reader.readScanline();
		
		typeConvert(scanlineIn, type, scanlineOut, UShort,
			dim.vec[0] * numChannels);

		// get scanline from the ImageMagick image for write access
#if defined (_WIN32) || defined (_WIN64) 

		Quantum* pixels = view.get(0, i, dim.vec[0], 1);

		for (j = 0; j < dim.vec[0]; j++, pixels++) {
			
			*pixels++ = scanlineOut[j*numChannels + channels.vec[0]];

			if (numOutChannels == 1 || numOutChannels == 2)
			{
				*pixels++ = scanlineOut[j*numChannels + channels.vec[0]];

				*pixels++ = scanlineOut[j*numChannels + channels.vec[0]];
			}
			else
			{
				*pixels++ = scanlineOut[j*numChannels + channels.vec[1]];
				if (numOutChannels == 4)
					*pixels++ = scanlineOut[j*numChannels + channels.vec[2]];
				else
					*pixels = scanlineOut[j*numChannels + channels.vec[2]];
			}

			if (numOutChannels == 2)
				*pixels = scanlineOut[j*numChannels + channels.vec[1]];

			else if (numOutChannels == 4)
				*pixels = scanlineOut[j*numChannels + channels.vec[3]];
			//else *pixels = 0; 

		}
		
#else 
		
		PixelPacket *pixels = view.get(0, i, dim.vec[0], 1);

		// now copy all the channels for all pixels in the scanline
		for (j = 0; j< dim.vec[0]; j++, pixels++)
		{
			pixels->red = scanlineOut[j*numChannels + channels.vec[0]];
			if (numOutChannels == 1 || numOutChannels == 2)
			{
				pixels->green = scanlineOut[j*numChannels + channels.vec[0]];
				pixels->blue = scanlineOut[j*numChannels + channels.vec[0]];
			}
			else
			{
				pixels->green = scanlineOut[j*numChannels + channels.vec[1]];
				pixels->blue = scanlineOut[j*numChannels + channels.vec[2]];
			}
			if (numOutChannels == 2)
				pixels->opacity = scanlineOut[j*numChannels + channels.vec[1]];
			else if (numOutChannels == 4)
				pixels->opacity = scanlineOut[j*numChannels + channels.vec[3]];
			else pixels->opacity = 0;
		}
#endif

		// update the actual image with the scanline we just copied
		view.sync();
	}

	reader.disconnect();
	outImage.write(argv[argc - 1]);

	return 0;
}
