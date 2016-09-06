
#if defined (_WIN32) || defined (_WIN64)
#include <Windows.h>
#include <stdlib.h>
#include <cctype>
#include <algorithm>
#else
#include <unistd.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <sstream>
#include <list>
#include <vector>

#if defined (_WIN32) || defined (_WIN64)

#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#if defined (_WIN32) || defined (_WIN64)
#define usleep(x) Sleep(x/1000)
#endif

#include "MDA/Base/CommandlineParser.hh"
#include "MDA/Expressions/ExpressionParseTree.hh"
#include "MDA/Expressions/ExpressionHelperFunctions.hh"
#include "MDA/Array/MDAFileIO.hh"


using namespace std;
using namespace MDA;


#define DEBUG_FLAG { cout << "DEBUG FLAG at " << __LINE__ << endl; }
#define USAGE "[options]"
#define TAB_WIDTH 4

#if defined (_WIN32) || defined (_WIN64)
#define MIN(x,y) x < y ? x: y 
#define MAX(x,y) x > y ? x: y 
#else
#define MIN(x, y) ({typeof(x) __x = (x); typeof(y) __y = (y); x < y ? x : y; })
#define MAX(x, y) ({typeof(x) __x = (x); typeof(y) __y = (y); x > y ? x : y; })
#endif

template <class T>
bool from_string(T& t, 
                 const string& s, 
                 ios_base& (*f)(std::ios_base&) = dec)
{
  istringstream iss(s);
  return !(iss >> f >> t).fail();
}


// Commands and options.
CommandlineParser parser;
list<CommandlineOption *> options;


// The input MDA.
MDAReader reader;
CoordinateVector dimensions;
DataType data_type;
unsigned int data_size;
unsigned int num_channels;
unsigned int num_scanlines;
unsigned char *raw_data; // I would like this to be void, but I can't do arith with that.

// Slices and animation

CoordinateVector slice_coord; // From and 3rd coord on.
bool playing;
long last_frame;
float fps = 24;

// Data manipulations.
EXPR::ExpressionSequence variable_exprs;
EXPR::ExpressionSequence channel_exprs;

// Positioning state.
int window_w = 400;
int window_h = 300;
int img_w;
int img_h;
float scale = 1;
float offset_x = 0;
float offset_y = 0;

// Mouse state.
int last_mouse_x;
int last_mouse_y;
int last_click_time = 0;
int mouse_modifiers;

#define FONT GLUT_BITMAP_8_BY_13
#define FONT_W 8
#define FONT_H 13

// Shell state.
bool sh_visible = false;
string sh_buffer; // The current command being typed.
int sh_carat_pos; // Location of insertion point.
list<string> sh_output; // All previous commands and output.
list<string> sh_history; // Commands we can jump to.
list<string>::iterator sh_history_it;
list<string> command_names;


bool show_coords = true;
bool show_values = true;
bool show_raw_data = false;
bool show_variables = false;


#define PS1 "$ "
#define SHELL_OUTPUT_LENGTH 100


// Background.
GLubyte bg_checker[128];
float bg_colours[] = {0.1, 0.25}; // {0.9, 0.75};

// What scales to start/finish fading in the pixel grid.
float grid_thresholds[] = {5, 20};


#define REGISTER_OPTION(name, T, var, msg, long, short, ...) { \
    bool has_short = strlen(short); \
    parser.registerOption(new T(var, msg, "--" long, has_short ? "-" short : NULL, ##__VA_ARGS__)); \
    options.push_back(new T(var, msg, long, has_short ? short : NULL, ##__VA_ARGS__)); \
    command_names.push_back(name); \
}


// Some prototypes.
void clean_fresh_options();
void animate();
void setup_animation_callback();


//
// some string processing ops
//

// strip whitespace from front and back of an STL string
std::string
strip( const std::string &str )
{
  unsigned long length= str.size();
  unsigned long i, j;
  
  for( i= 0 ; i< length ; i++ )
    if( !std::isspace( str[i] ) )
      break;
  
  for( j= length-1 ; j> i ; j-- )
    if( !std::isspace( str[j] ) )
      break;
  
  if( i== 0 && j== length-1 )
    return str;
  else
    return str.substr( i, (j+1)-i );
}


// replace all instances of a substring with a diffeent string
std::string
replace( const std::string &str,
	 const std::string &oldS, const std::string &newS )
{
  unsigned long oldLen= oldS.size();
  long oldPos= 0;
  std::string result( str );
  
  while( (oldPos= result.find( oldS, oldPos ))>= 0 )
  {
    cerr << oldPos << endl;
    result.replace( oldPos, oldLen, newS );
  }
  
  return result;
}

// split string at a provided separator character
void
split( const std::string &str, std::vector<std::string> &result, char sep )
{
  unsigned long length= str.size();
  unsigned long i, j;
  
  for( i= j= 0 ; j< length ; j++ )
    if( str[j]== sep )
    {
      result.push_back( str.substr( i, j-i ) );
      i= j+1;
    }
  result.push_back( str.substr( i, length-i ) );
}
  

void init_main(int argc, char *argv[]) {
    
    // Build the background checker stipple.
    for (int i = 0; i < 32; i++) {
        bg_checker[4 * i + 0] = i < 16 ? 0 : 0xFF;
        bg_checker[4 * i + 1] = i < 16 ? 0 : 0xFF;
        bg_checker[4 * i + 2] = i < 16 ? 0xFF : 0;
        bg_checker[4 * i + 3] = i < 16 ? 0xFF : 0;
    }
    
    REGISTER_OPTION("variables, v", ExpressionOption, variable_exprs,
         "\tdefinitions of temporary variables (a sequence of expressions)\n",
         "variables", "v"
     );
    
    REGISTER_OPTION("channels, c", ExpressionOption, channel_exprs,
         "\tdefinitions of temporary variables (a sequence of expressions)\n",
         "channels", "c"
     );
    
    REGISTER_OPTION("scale, s", FloatOption, scale,
        "\tscaling factor of the image\n",
        "scale", "s",
        0.01f, 100.0f
    );
    
    REGISTER_OPTION("slice", CoordinateOption, slice_coord,
        "\tslice of the image stack to display\n",
        "slice", ""
    )
    
    REGISTER_OPTION("play, stop", BoolOption, playing,
        "\tplay stack as an animation\n",
        "play", "",
        "stop", ""
    );
    
    REGISTER_OPTION("fps", FloatOption, fps,
        "\tframe rate of stack animation\n",
        "fps", "",
        0.01f, 100.0f
    );
    
    REGISTER_OPTION("show-coords, hide-coords", BoolOption, show_coords,
        "\tdisplay coordinates in overlay\n",
        "show-coords", "",
        "hide-coords", ""
    );
    
    REGISTER_OPTION("show-values, hide-values", BoolOption, show_values,
        "\tdisplay pixel values in overlay\n",
        "show-values", "",
        "hide-values", ""
    );
    
    REGISTER_OPTION("show-raw-data, hide-raw-data", BoolOption, show_raw_data,
        "\tdisplay raw file data in overlay\n",
        "show-raw-data", "",
        "hide-raw-data", ""
    );
    
    REGISTER_OPTION("show-variables, hide-variables", BoolOption, show_variables,
        "\tdisplay variables in overlay\n",
        "show-variables", "",
        "hide-variables", ""
    );
    
    // parse options
    int index = 1;
    if(!parser.parse(index, argc, argv)) {
        parser.usage(argv[0], USAGE);
        exit(1);
    }
    
    // Open the MDA via stdin.
    reader.connect(cin);
    if (!reader.readHeader()) {
        cerr << "Cannot read file header." << endl;
        exit(3);
    }
    data_type = reader.getType();
    data_size = dataTypeSizes[data_type];
    dimensions = reader.getDim();
    num_channels = reader.getNumChannels();
    num_scanlines = reader.getNumScanlinesLeft();
    
    // Needs to be after we load the MDA so we know how many channels we want.
    clean_fresh_options();
                
    if (dimensions.vec.size() < 2) {
        cerr << "MDA must have atleast 2 dimensions." << endl;
        exit(4);
    }
    
    img_w = dimensions.vec[0];
    img_h = dimensions.vec[1];
    
    // Alloc a pointer for all of the image data.
    raw_data = (unsigned char*)malloc(num_scanlines * img_w * num_channels * data_size);
    
    // Read in full image data and convert to floats.
    // Maybe should be doing this to double... oops.
    for (int i = 0; i < num_scanlines; i++) {
        void *scanline = reader.readScanline();
        memcpy(raw_data + (i * img_w * num_channels * data_size), scanline, img_w * num_channels * data_size);
    }
    
    // Prime the frame rate limiter;
    last_frame = glutGet(GLUT_ELAPSED_TIME);
}


void clean_fresh_options() {
    if (channel_exprs.size() == 0) {
        stringstream ss (stringstream::in | stringstream::out);
        for (int i = 0; i < min(4, (int)num_channels); i++) {
            ss << (i ? "," : "");
            ss << "#";
            ss << i;
        }
        channel_exprs = EXPR::parse(const_cast<char*>(ss.str().c_str()));
    }
    if (channel_exprs.size() > 4) {
        cerr << "Too many channels; truncating to 4." << endl;
        channel_exprs.erase(channel_exprs.begin() + 4);
    }
    
    // Assert that the slice coord is of the appropriate dimension.
    slice_coord.vec.resize(dimensions.vec.size() - 2, 0);
    
    setup_animation_callback();
}


// Returns offset into raw data (ignoring data type) given the x/y coords, and
// the remaining coords from slice_coord.
static inline int offset_into_raw_data(int x, int y) {
    int dim_size = num_channels;
    int offset = 0;
    int coord;
    for (int i = 0; i < dimensions.vec.size(); i++) {
        switch (i) {
            case 0: coord = x; break;
            case 1: coord = y; break;
            default: coord = slice_coord.vec[i - 2]; break;
        }
        offset += dim_size * coord;
        dim_size *= dimensions.vec[i];
    }
    return offset;
}


static void get_raw_pixel_values(int x, int y, double *out) {
    typeConvert(raw_data + offset_into_raw_data(x, y) * data_size, data_type, (void*)out, Double, num_channels); 
}


static void get_expr_pixel_values(int x, int y, float *val_out, float *var_out = NULL) {
#if defined (_WIN32) || defined (_WIN64)
	double* raw_data = new double [num_channels];
	double* variables = new double [variable_exprs.size()];
#else

    double raw_data[num_channels];
    double variables[variable_exprs.size()];
#endif
    get_raw_pixel_values(x, y, raw_data);
    for (int k = 0; k < variable_exprs.size(); k++) {
        variables[k] = variable_exprs[k]->eval(raw_data, variables);
        if (var_out) {
            var_out[k] = variables[k];
        }
    }
    for (int k = 0; k < channel_exprs.size(); k++) {
        val_out[k] = channel_exprs[k]->eval(raw_data, variables);
    }
#if defined (_WIN32) || defined (_WIN64)
	delete[] raw_data;
	delete[] variables;
#endif

}


static void init_opengl() {
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}


static void init_texture() {

    float *tex_data = new float[img_w * img_h * channel_exprs.size()];

    // Pixel by pixel vxaluate all pixels to build the texture.
    for (int x = 0; x < img_w; x++) {
        for (int y = 0; y < img_h; y++) {
            // Not offset_into_raw_data.
            int offset = (y * img_w + x) * channel_exprs.size();
            get_expr_pixel_values(x, y, tex_data + offset);
        }
    }
    
    // Bind texture
    static int channel_types[] = {0, GL_LUMINANCE, GL_LUMINANCE_ALPHA, GL_RGB, GL_RGBA};
    int channel_type = channel_types[channel_exprs.size()];
    glTexImage2D(
        GL_TEXTURE_2D,
        0,
        channel_type,
        img_w,
        img_h,
        0,
        channel_type,
        GL_FLOAT,
        tex_data
    );

    delete [] tex_data;
}


static void reshape(int w, int h) {
    window_w = w;
    window_h = h;
    
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, w, 0, h, -1, 1);
    glMatrixMode(GL_MODELVIEW);

}



static void zoom(float amount) {
    scale *= amount;
    offset_x *= amount;
    offset_y *= amount;
}



static void add_to_shell_output(string input) {
    size_t start = 0, pos = 0;
    do {
        pos = input.find('\n', start);
        if (pos == string::npos) {
            sh_output.push_back(input.substr(start));
        } else {
            sh_output.push_back(input.substr(start, pos-start));
            start = pos + 1;
        }
    } while(pos != string::npos);
}





static void shell_exec(const string &line) {
    
    int argc;
    char *argv[2];
    string command, argv_command, arg;
    list<CommandlineOption *>::iterator opt_it;
    
    int i = line.find(' ');
    if (i == string::npos) {
        argc = 1;
        command = line;
    } else {
        argc = 2;
        command = line.substr(0, i);
        arg = line.substr(i + 1);
    }
    argv[0] = const_cast<char*>(command.c_str());
    argv[1] = const_cast<char*>(arg.c_str());
        
    if (command == "help") {        
        if (argc == 1) {
            stringstream ss;
            ss << "availible commands (try \"help <command>\"):\n";
            list<string>::iterator name_it;
            for (name_it = command_names.begin(); name_it != command_names.end(); name_it++) {
                ss << "\t" << *name_it << "\n";
            }
            add_to_shell_output(ss.str());
        } else {
            bool found = false;
            for (opt_it = options.begin(); opt_it != options.end(); opt_it++) {
                if ((*opt_it)->isResponsible(argv[1])) {
                    stringstream ss;
                    (*opt_it)->usage(ss);
                    add_to_shell_output(ss.str());
                    found = true;
                    break;
                }
            }
            if (!found) {
                add_to_shell_output("command not found: \"" + arg + "\"");
            }
        }
    }
    
    else if ((command == "c" || command == "channels") && arg.find("$") != string::npos) {
        
        for (int i = 0; i < channel_exprs.size(); i++) {
            stringstream ss;
            ss << "(";
            channel_exprs[i]->printTree(ss);
            ss << ")";
            string new_expr = replace(arg, "$", ss.str());
            delete channel_exprs[i];
            // XXX: memory leak on more than one new expr.
            channel_exprs[i] = EXPR::parse(const_cast<char*>(new_expr.c_str()))[0];
        }
        init_texture();
        
        
    }
    /* // Replace by channel overloading; will be removed soon
    else if (command == "pow") {
        if (argc == 1) {
            add_to_shell_output("usage: pow <exponent>");
            return;
        }
        float factor;
        if (!from_string<float>(factor, arg)) {
            add_to_shell_output("\"" + arg + "\" is not a float");
            return;
        }
        for (int i = 0; i < channel_exprs.size(); i++) {
            ExpressionParseTree *factor_expr = new ExpressionParseTree;
            ExpressionParseTree *func_expr = new ExpressionParseTree;
            factor_expr->makeConstant(factor);
            func_expr->makeBivariate(&pow, channel_exprs[i], factor_expr);
            delete channel_exprs[i];
            channel_exprs[i] = func_expr;
        }
        init_texture();
    }
    
    else if (command == "invert") {
        for (int i = 0; i < channel_exprs.size(); i++) {
            ExpressionParseTree *one_expr = new ExpressionParseTree;
            ExpressionParseTree *sub_expr = new ExpressionParseTree;
            one_expr->makeConstant(1);
            sub_expr->makeBivariate(&minusOp, one_expr, channel_exprs[i]);
            delete channel_exprs[i];
            channel_exprs[i] = sub_expr;
        }
        init_texture();
    }
    //*/
    
    else {
        
    
        // Check all of the CommandlineOptions to see if we have a match.
        int index = 1;
        bool found = false;
        for (opt_it = options.begin(); opt_it != options.end(); opt_it++) {
            if ((*opt_it)->isResponsible(argv[0])) {
                if (!(*opt_it)->parse(index, argc, argv)) {
                    stringstream ss;
                    (*opt_it)->usage(ss);
                    add_to_shell_output(ss.str());
                }
                found = true;
                break;
            }
        }
        if (found) {
            clean_fresh_options();
            init_texture();
        } else {
            add_to_shell_output("command not found: \"" + command + "\"; try \"help\"");
        }
    }
}

static void special_keyboard(int key, int x, int y);
static void keyboard(unsigned char key, int x, int y) {
    
    glutPostRedisplay();
    
    if (key == 27) { // escape
        sh_visible = !sh_visible;
        return;
    }
    
    if (sh_visible) {
        
        // Mutate keys into other keys.
        if (key == 9) { // Tab
            key = ' ';
        }
        
        // Variables for shell execution
        vector<string> sh_commands;
        
        // Handle all control characters:        
        switch(key) {
        
        case 13: // return
        case 10: // ctrl-J
            
            sh_buffer = strip(sh_buffer);
            sh_output.push_back(PS1 + sh_buffer);
            if (sh_history_it != sh_history.end() && !sh_history.empty()) {
                sh_history.pop_back();
            }
            sh_history.push_back(sh_buffer);
            
            // Split on semicolons.
            split(sh_buffer, sh_commands, ';');
            for (int i = 0; i < sh_commands.size(); i++) {
                string sh_command = strip(sh_commands[i]);
                if (sh_command.size()) {
                    shell_exec(sh_command);
                }
            }

            // Don't want the output/history getting too big.
            while (sh_output.size() > SHELL_OUTPUT_LENGTH) {
                sh_output.pop_front();
            }
            while (sh_history.size() > SHELL_OUTPUT_LENGTH) {
                sh_history.pop_front();
            }
            
            sh_buffer = "";
            sh_carat_pos = 0;
            sh_history_it = sh_history.end();
            break;
            
        case 127: // backspace
            sh_buffer = sh_buffer.substr(0, sh_carat_pos - 1) + sh_buffer.substr(sh_carat_pos);
            sh_carat_pos = MAX(0, sh_carat_pos - 1);
            break;
            
        case 4:
        case 8: // delete | ctrl-H
            if (sh_carat_pos != sh_buffer.size()) {
                sh_buffer = sh_buffer.substr(0, sh_carat_pos) + sh_buffer.substr(sh_carat_pos + 1);
                sh_carat_pos = MIN(sh_carat_pos, sh_buffer.size());
            }
            break;
        
        case 1: // ctrl-A
            sh_carat_pos = 0;
            break;
        
        case 5: // ctrl-E
            sh_carat_pos = sh_buffer.size();
            break;
        
        case 6: // ctrl-F
            special_keyboard(GLUT_KEY_RIGHT, x, y);
            break;
            
        case 2: // ctrl-B
            special_keyboard(GLUT_KEY_LEFT, x, y);
            break;
            
        case 16: // ctrl-P
            special_keyboard(GLUT_KEY_UP, x, y);
            break;
            
        case 14: // ctrl-N
            special_keyboard(GLUT_KEY_DOWN, x, y);
            break;
            
        case 11: // ctrl-K
            if (sh_carat_pos != sh_buffer.size()) {
                sh_buffer = sh_buffer.substr(0, sh_carat_pos);
            }
            break;
        
        case 21: // ctrl-U
            sh_buffer = "";
            sh_carat_pos = 0;
            break;
        
        default:
            // unknown control characters
            if (key < 32) {
                cerr << "unknown control character: " << int(key) << endl;
            }
            // printable characters
            else {
                sh_buffer = sh_buffer.substr(0, sh_carat_pos) + string(1, key) + sh_buffer.substr(sh_carat_pos);
                sh_carat_pos += 1;
            }
        }
        
    } else {
        
        switch(key) {
    
        case ' ':
            playing = !playing;
            setup_animation_callback();
            break;
        case 'a':
            zoom(2);
            break;
        case 'z':
            zoom(0.5);
            break;
        case 'r':
            scale = 1;
            offset_x = 0;
            offset_y = 0;
            break;
        }
        
    }
}


static void special_keyboard(int key, int x, int y) {
    glutPostRedisplay();
    if (sh_visible) {
        switch (key) {
            
            // Moving the carat.
            case GLUT_KEY_LEFT:
                sh_carat_pos = MAX(0, sh_carat_pos - 1);
                break;
            case GLUT_KEY_RIGHT:
                sh_carat_pos = MIN(sh_carat_pos + 1, sh_buffer.size());
                break;
                
            // Moving through history.
            case GLUT_KEY_UP:
                if (sh_history_it != sh_history.begin() && !sh_history.empty()) {
                    if (sh_history_it == sh_history.end() && sh_history.size()) {
                        // We must do this as well otherwise we will shift
                        // with the end of the list.
                        --sh_history_it;
                        sh_history.push_back(sh_buffer);
                    } else {
                    --sh_history_it;
                    }
                    sh_buffer = *sh_history_it;
                    sh_carat_pos = sh_buffer.size();
                }
                break;
            case GLUT_KEY_DOWN:
                if (sh_history_it != sh_history.end()) {
                    ++sh_history_it;
                    if (sh_history_it == sh_history.end()) {
                        sh_buffer = "";
                    } else {
                        sh_buffer = *sh_history_it;
                    }
                    sh_carat_pos = sh_buffer.size();
                }
                break;
            
        }
    } else {
        int slice_i;
        int direction = 0;
        switch(key) {
            case GLUT_KEY_RIGHT:
                direction = 1;
                slice_i = 0;
                break;
            case GLUT_KEY_LEFT:
                direction = -1;
                slice_i = 0;
                break;
            case GLUT_KEY_UP:
                direction = 1;
                slice_i = 1;
                break;
            case GLUT_KEY_DOWN:
                direction = -1;
                slice_i = 1;
                break;
        }
        if (direction && slice_i < slice_coord.vec.size()) {
            // A little extra protection here because negative mods are
            // implementation specific.
            slice_coord.vec[slice_i] = (slice_coord.vec[slice_i] + dimensions.vec[slice_i + 2] + direction) % dimensions.vec[slice_i + 2];
            init_texture();
        }
    }
}


static void mouse_to_img(double x, double y, double *img_x, double *img_y) {
    double model_mat[16];
    double proj_mat[16];
    int viewport[4];
    double out[3];
    glGetDoublev(GL_PROJECTION_MATRIX, proj_mat);
    glGetDoublev(GL_MODELVIEW_MATRIX, model_mat);
    glGetIntegerv(GL_VIEWPORT, viewport);
    gluUnProject(
        x, window_h - y, 0, // Need to invert the y coords to match OpenGL.
        model_mat, proj_mat, viewport,
        out, out + 1, out + 2
    );
    *img_x = out[0];
    *img_y = out[1];
}


static void mouse_double(int, int, int);

static void mouse(int button, int state, int x, int y) {
    
    // button is one of GUT_(LEFT|MIDDLE|RIGHT)_BUTTON
    // state is one of GLUT_(UP|DOWN)
    // glutGetModifiers() will return ORed GLUT_ACTIVE_(SHIFT|CTRL|ALT)
        
    if (state == GLUT_DOWN) {
        mouse_modifiers = glutGetModifiers();
        
        int click_time = glutGet(GLUT_ELAPSED_TIME);
        if (click_time - last_click_time < 250) {
            mouse_double(button, x, y);
            last_click_time = 0; // Don't want triple clicks to count.
        }
        else {
            last_click_time = click_time;
        }
    }
    else {
        mouse_modifiers = 0;
    }
    
    last_mouse_x = x;
    last_mouse_y = y;
}


static void mouse_double(int button, int x, int y) {
    
    // Calculate the distance from the center of the image to the click.
    int dx = x - (window_w / 2 + offset_x);
    int dy = y - (window_h / 2 + offset_y);
    
    // We want to manipulate scale and offsets so that dx will be twice what
    // it was before. Derivation:
    // dx = x - (w/2 + o) -> o = x - dx - w/2
    // set dx' = 2dx:
    // o' = x - 2dx - w/2
    offset_x = x - 2 * dx - window_w / 2;
    offset_y = y - 2 * dy - window_h / 2;
    scale *= 2;
    
    glutPostRedisplay();
    
}


static void mouse_active_move(int x, int y) {
    
    if (mouse_modifiers & GLUT_ACTIVE_ALT) {
        int dx = x - last_mouse_x;
        float factor = 1 + ((float)dx / 150.0f);
        zoom(factor);
    }
    else {
        offset_x += x - last_mouse_x;
        offset_y += y - last_mouse_y;
    }
    
    glutPostRedisplay();
    last_mouse_x = x;
    last_mouse_y = y;
}


static void mouse_passive_move(int x, int y) {
    glutPostRedisplay();
    last_mouse_x = x;
    last_mouse_y = y;
}






template <class T>
static void print(T x, int line_margin = 0) {
    stringstream ss;
    ss << x;
    print(ss.str());
}

template <>
void print <char> (char c, int line_margin) {
    if (c == '\t') {
        for (int i = 0; i < TAB_WIDTH; i++) {
            glutBitmapCharacter(FONT, ' ');
        }
    } else {
        glutBitmapCharacter(FONT, c);
    }
}

template <>
void print <string> (string s, int line_margin) {
    int size = s.size();
    for (int i = 0; i < size; i++) {
        print(s[i]);
    }
}




static void display(void) {
    
    glLoadIdentity();
    
    // Draw the background.
    glClearColor(bg_colours[0], bg_colours[0], bg_colours[0], 1);
    glClear(GL_COLOR_BUFFER_BIT);
    glPolygonStipple(bg_checker);
    glEnable(GL_POLYGON_STIPPLE);
    glColor3f(bg_colours[1], bg_colours[1], bg_colours[1]);
    glBegin(GL_POLYGON);
    glVertex2f(0, 0);
    glVertex2f(window_w, 0);
    glVertex2f(window_w, window_h);
    glVertex2f(0, window_h);
    glEnd();
    glDisable(GL_POLYGON_STIPPLE);
    
    // Position the image.
    glTranslatef(window_w / 2, window_h / 2, 0); // This should be int division.
    glScalef(1, -1, 1); // To match GLUT/image coordinates.
    glTranslatef(offset_x, offset_y, 0);
    glScalef(scale, scale, 1);
    glTranslatef(-img_w / 2, -img_h / 2, 0);

    // Draw the image.
    glEnable(GL_TEXTURE_2D);
    glColor3f(1, 1, 1);
    glBegin(GL_POLYGON);
    glTexCoord2f(0, 0); glVertex2f(0, 0);
    glTexCoord2f(1, 0); glVertex2f(img_w, 0);
    glTexCoord2f(1, 1); glVertex2f(img_w, img_h);
    glTexCoord2f(0, 1); glVertex2f(0, img_h);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    
    // Get mouse/img coords.
    double mouse_img_x, mouse_img_y;
    mouse_to_img(last_mouse_x, last_mouse_y, &mouse_img_x, &mouse_img_y);
    int mx = mouse_img_x >= 0 ? mouse_img_x : -1;
    int my = mouse_img_y >= 0 ? mouse_img_y : -1;
    
    // Overlay pixel values.
    if (mx >= 0 && mx < img_w && my >= 0 && my < img_h) {
        
        glPushMatrix();
        glLoadIdentity();
        glTranslatef(10, window_h - 5, 0);
        
        #define NEWLINE(margin) {glTranslatef(0, -FONT_H - margin, 0); glRasterPos2f(0, 0);}
        NEWLINE(0);
        
        // Coordinates.
        if (show_coords) {
            print(' '); print(mx); print(','); print(my);
            if (slice_coord.vec.size()) {
                print(',');
                print(slice_coord);
            }
            NEWLINE(5);
        }
        
        // Final colours.
        static const char *channel_names[] = {"L", "LA", "RGB", "RGBA"};
		#if defined (_WIN32) || defined (_WIN64)
		float* values=new float[channel_exprs.size()];
		float* variables=new float[variable_exprs.size()];
#else
        float values[channel_exprs.size()];
        float variables[variable_exprs.size()];
#endif
        get_expr_pixel_values(mx, my, values, variables);
        if (show_values) {
            for (int i = 0; i < channel_exprs.size(); i++) {
                print(' ');
                print(channel_names[channel_exprs.size() - 1][i]);
                print(": ");
                print(values[i]);
                NEWLINE(0);
            }
            NEWLINE(-5);
        }
        
        // Variables
        if (show_variables) {
            for (int i = 0; i < variable_exprs.size(); i++) {
                print('%');
                print(i);
                print(": ");
                print(variables[i]);
                NEWLINE(0);
            }
            if (variable_exprs.size()) {
                NEWLINE(-5);
            }
#if defined (_WIN32) || defined (_WIN64)
			delete[] values;
			delete[] variables;
#endif
        }
        
        // Raw values.
        if (show_raw_data) {
#if defined (_WIN32) || defined (_WIN64)
			double* raw_values = new double [num_channels];
#else
            double raw_values[num_channels];
#endif 
            get_raw_pixel_values(mx, my, raw_values);
            for (int i = 0; i < num_channels; i++) {
                print('#');
                print(i);
                print(": ");
                print(raw_values[i]);
                NEWLINE(0);
            }

#if defined (_WIN32) || defined (_WIN64)
			delete[] raw_values;
#endif
        }
        
        #undef NEWLINE
        glPopMatrix();
    }
    
    float grid_opac = (scale - grid_thresholds[0]) / (grid_thresholds[1] - grid_thresholds[0]);
    if (grid_opac > 0) {
        grid_opac = min(1.0f, grid_opac);
        
        // Pixel highlight.
        if (mx >= 0 && mx < img_w && my >= 0 && my < img_h) {
            glColor4f(1, 1, 0, 0.25 * grid_opac);
            glBegin(GL_POLYGON);
            glVertex2f(mx, my);
            glVertex2f(mx + 1, my);
            glVertex2f(mx + 1, my + 1);
            glVertex2f(mx, my + 1);
            glEnd();
        }
        
        // Grid lines.
        glColor4f(1, 1, 1, grid_opac);
        for (int x = 0; x <= img_w; x++) {
            glBegin(GL_LINES);
            glVertex2f(x, 0);
            glVertex2f(x, img_h);
            glEnd();
        }
        for (int y = 0; y <= img_h; y++) {
            glBegin(GL_LINES);
            glVertex2f(0, y);
            glVertex2f(img_w, y);
            glEnd();
        }
        
    }
    
    // Shell
    if (sh_visible) {
        glPushMatrix();
        glLoadIdentity();
        
        // Background.
        glColor4f(0, 0, 0, .65);
        glBegin(GL_POLYGON);
        glVertex2f(0, 0);
        glVertex2f(0, window_h);
        glVertex2f(window_w, window_h);
        glVertex2f(window_w, 0);
        glEnd();
        
        glTranslatef(5, 5, 0);
        
        #define NEWLINE(margin) {glTranslatef(0, FONT_H + margin, 0); glRasterPos2f(0, 0);}
        
        // Insertion carat.
        glPushMatrix();
        glTranslatef((2 + sh_carat_pos) * FONT_W, 0, 0);
        glColor4f(1, 1, 1, 0.25);
        glBegin(GL_POLYGON);
        glVertex2f(0, 0); glVertex2f(0, FONT_H); glVertex2f(FONT_W, FONT_H); glVertex2f(FONT_W, 0);
        glEnd();
        glPopMatrix();
        
        glColor4f(1, 1, 1, 0.75);
        glRasterPos2f(0, 0);
        
        // Print the buffer.
        // TODO: line wrapping.
        string buffer = PS1 + sh_buffer;
        print(buffer);
        NEWLINE(0);
        
        // Print the output.
        // TODO: line wrapping.
        list<string>::reverse_iterator rit;
        for (rit = sh_output.rbegin(); rit != sh_output.rend(); ++rit) {
            print(*rit);
            NEWLINE(0);
        }
        
        #undef NEWLINE
        glPopMatrix();
    }
    
    glutSwapBuffers();
}


void setup_animation_callback() {
    glutIdleFunc(playing && slice_coord.vec.size() ? animate : NULL);
}

void animate(void) {
    
    // Rate limiting to the frame rate.
    int current_time = glutGet(GLUT_ELAPSED_TIME);
    int next_frame = last_frame + (1000.0f / fps);
    if (current_time < next_frame) {
        usleep(MIN(5000, 1000 * (next_frame - current_time)));
        return;
    }
    last_frame = next_frame;
    
    // Cycle the slice coord
    slice_coord.vec[0] = (slice_coord.vec[0] + 1) % dimensions.vec[2];
    init_texture();
    glutPostRedisplay();
}


int main(int argc, char *argv[]) {
    
    glutInit(&argc, argv);
    init_main(argc, argv);
        
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(min(1024, img_w), min(768, img_h));
    
    glutCreateWindow(argv[0]);
    
    init_opengl();
    init_texture();
    
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special_keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(mouse_active_move);
    glutPassiveMotionFunc(mouse_passive_move);
    
    glutMainLoop();
    return 0;
    
}
