/**
 * Reads and shows 3D MDA files
 * 
 * @author Mike Krimerman <krim@cs.ubc.ca>
 */

import java.io.File;
import java.io.RandomAccessFile;
import java.io.ByteArrayInputStream;
import java.io.DataInputStream;

import ij.io.OpenDialog;
import ij.plugin.*;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;
import ij.IJ;

public class Open_MDA extends ImagePlus implements PlugIn 
{
    public void run( String path )
    {
	boolean needToShow = false;
	File f;

	if ( path == null || path.equals("") ) 
	{
	    OpenDialog dialog = new OpenDialog( "Open a volumetric MDA file", null );
	    if ( dialog.getDirectory() == null )
		return; // canceled
	    f = new File( dialog.getDirectory(), dialog.getFileName() );
	    /*
	     * Since no path was passed, assume that it was run interactively
	     * rather than from HandleExtraFileTypes
	     */
	    needToShow = true;
	}
	else 
	{
	    f = new File( path );
	}

	try 
	{
/*
MDA 1.0
Format:         float little
Dimensions:     1920,1080
Channels:       8
###
*/
	    RandomAccessFile file = new RandomAccessFile( f, "r" );
	    String version = file.readLine();	
	    if ( !version.equals( "MDA 1.0" ) )	
	    {
		IJ.error( "'" + f + "' is not an MDA file\n" );
		return;
	    }
	    String[] tokens = file.readLine().split( "[\t ]+" );
	    String mdaType = tokens[ 1 ];
	    String mdaEndian = tokens[ 2 ];

	    tokens = file.readLine().split( "[\t ]+" ); 
	    String[] mdaDimensions = tokens[ 1 ].split( "," );
	    
	    tokens = file.readLine().split( "[\t ]+" ); 
	    String mdaChannels = tokens[ 1 ];

	    String terminator = file.readLine(); // Hash marks, end of header

	    boolean littleEndian = mdaEndian.equals( "little" );
	    int numChannels = Integer.parseInt( mdaChannels );
	    if ( numChannels != 1 ) 
	    {
		IJ.error( "Unsupported number of channels (" + Integer.toString( numChannels ) + "), must contain a single channel\n" );
		return;
	    }
	    if ( mdaDimensions.length < 2 || mdaDimensions.length > 3 ) 
	    {
		IJ.error( "Unsupported number of dimensions (" + Integer.toString( mdaDimensions.length ) + "), must be 2 or 3\n" );
		return;
	    }
	    int width  = Integer.parseInt( mdaDimensions[ 0 ] );
	    int height = Integer.parseInt( mdaDimensions[ 1 ] );
	    int depth  = mdaDimensions.length == 3 ? Integer.parseInt( mdaDimensions[ 2 ] ) : 1;

	    byte typeLength = 0;
	    if ( mdaType.equals( "ubyte" ) ) 
	    {
		typeLength = 1;
	    }
	    else if ( mdaType.equals( "ushort" ) ) 
	    {
		typeLength = 2;
	    }
	    else if ( mdaType.equals( "float" ) ) 
	    {
		typeLength = 4;
	    }
			
	    if ( typeLength == 0 )
	    {
		IJ.error( "Unsupported MDA type (" + mdaType + "). Supported types: ubyte, ushort, float.\n" );
		return;
	    }
			
	    ImageStack stack = new ImageStack( width, height );
	    float min = Float.MAX_VALUE;
	    float max = Float.MIN_VALUE;
			
	    //IJ.log( "MDA " + width + "x" + height + "x" + depth + ", " + mdaType + " " + mdaEndian );
	    for ( int z = 0; z < depth; ++z )
	    {
		ImageProcessor ip = null;
		
		switch ( typeLength )
		{
		case 1:
		    byte[] byteSlice = new byte[ height * width ];
		    file.read( byteSlice );
		    ip = new ByteProcessor( width, height, byteSlice, null );
		    break;
		case 2:
		    {
			byte[] byteBuffer = new byte[ height * width * 2 ];
			file.read( byteBuffer );
			if ( littleEndian )
			{
			    for ( int i = 0; i < height * width; ++i )
			    {
				byte b0 = byteBuffer[ 2*i + 0 ]; 
				byte b1 = byteBuffer[ 2*i + 1 ]; 
				byteBuffer[ 2*i + 0 ] = b1; 
				byteBuffer[ 2*i + 1 ] = b0; 
			    }
			}
			ByteArrayInputStream bais = new ByteArrayInputStream( byteBuffer );
			DataInputStream dis = new DataInputStream( bais );
			short[] shortSlice = new short[ height * width ];
			for ( int i = 0; i < shortSlice.length; ++i )
			{
			    short value = dis.readShort(); 
			    shortSlice[ i ] = value;
			    if ( value > max ) max = value;
			    if ( value < min ) min = value;
			}
			ip = new ShortProcessor( width, height, shortSlice, null );
		    }
		    break;
		case 4: 
		    {
			byte[] byteBuffer = new byte[ height * width * 4 ];
			file.read( byteBuffer );
			if ( littleEndian )
			{
			    for ( int i = 0; i < height * width; ++i )
			    {
				byte b0 = byteBuffer[ 4*i + 0 ]; 
				byte b1 = byteBuffer[ 4*i + 1 ]; 
				byte b2 = byteBuffer[ 4*i + 2 ]; 
				byte b3 = byteBuffer[ 4*i + 3 ]; 
				byteBuffer[ 4*i + 0 ] = b3; 
				byteBuffer[ 4*i + 1 ] = b2; 
				byteBuffer[ 4*i + 2 ] = b1; 
				byteBuffer[ 4*i + 3 ] = b0; 
			    }
			}
			ByteArrayInputStream bais = new ByteArrayInputStream( byteBuffer );
			DataInputStream dis = new DataInputStream( bais );
			float[] floatSlice = new float[ height * width ];
			for ( int i = 0; i < height * width; ++i )
			{
			    float value = dis.readFloat();
			    floatSlice[ i ] = value;
			    if ( value > max ) max = value;
			    if ( value < min ) min = value;
			}
			ip = new FloatProcessor( width, height, floatSlice, null );
		    }
		    break;
		}
		IJ.showProgress( z + 1, depth );
		IJ.showStatus( "Reading: " + (z + 1) + "/" + depth );
		stack.addSlice( null, ip );
	    }
	    ImagePlus imp = new ImagePlus( f.getName().replaceAll( ".mda$", "" ), stack );
	    imp.setDisplayRange( min, max );
	    
	    if ( needToShow )
	    {
		imp.show();
	    }
	}
	catch ( Exception e )
	{
	    IJ.error( "Opening '" + f + "' as mda failed.\n" + e.getMessage() );
	}
    }
}
