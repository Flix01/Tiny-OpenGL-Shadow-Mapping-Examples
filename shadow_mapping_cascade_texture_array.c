// https://github.com/Flix01/Tiny-OpenGL-Shadow-Mapping-Examples
/** License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/

// DEPENDENCIES:
/*
-> glut or freeglut (the latter is recommended)
-> glew (Windows only)
*/

// HOW TO COMPILE:
/*
// LINUX:
gcc -O2 -std=gnu89 shadow_mapping_cascade_texture_array.c -o shadow_mapping_cascade_texture_array -I"../" -lglut -lGL -lX11 -lm
// WINDOWS (here we use the static version of glew, and glut32.lib, that can be replaced by freeglut.lib):
cl /O2 /MT /Tc shadow_mapping_cascade_texture_array.c /D"GLEW_STATIC" /I"../" /link /out:shadow_mapping_cascade_texture_array.exe glut32.lib glew32s.lib opengl32.lib gdi32.lib Shell32.lib comdlg32.lib user32.lib kernel32.lib


// IN ADDITION:
By default the source file assumes that every OpenGL-related header is in "GL/".
But you can define in the command-line the correct paths you use in your system
for glut.h, glew.h, etc. with something like:
-DGLUT_PATH=\"Glut/glut.h\"
-DGLEW_PATH=\"Glew/glew.h\"
(this syntax works on Linux, don't know about Windows)
*/

//#define USE_GLEW  // By default it's only defined for Windows builds (but can be defined in Linux/Mac builds too)


#define PROGRAM_NAME "shadow_mapping_cascade_texture_array"
#define VISUALIZE_DEPTH_TEXTURE
//#define VISUALIZE_CASCADE_SPLITS
#define SHADOW_MAP_RESOLUTION 512
#define SHADOW_MAP_NUM_CASCADES 4
#define SHADOW_MAP_CASCADE_LAMBDA   0.7     // in [0=uniform splits,1=logarithmic splits] logarithmic splits put higher resolution near the camera
#define SHADOW_MAP_CLAMP_MODE GL_CLAMP_TO_EDGE // GL_CLAMP or GL_CLAMP_TO_EDGE or GL_CLAMP_TO_BORDER
    //          GL_CLAMP;               // sampling outside of the shadow map gives always shadowed pixels
    //          GL_CLAMP_TO_EDGE;       // sampling outside of the shadow map can give shadowed or unshadowed pixels (it depends on the edge of the shadow map)
    //          GL_CLAMP_TO_BORDER;     // sampling outside of the shadow map gives always non-shadowed pixels (if we set the border color correctly)
#define SHADOW_MAP_FILTER GL_LINEAR // GL_LINEAR or GL_NEAREST (GL_LINEAR is more useful with a sampler2DShadow, that cannot be used with esponential shadow mapping)

//#define USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE   // Better resolution, but shadow-swimming as the camera rotates (on static objects). Please see README.md about it.


// These path definitions can be passed to the compiler command-line
#ifndef GLUT_PATH
#   define GLUT_PATH "GL/glut.h"    // Mandatory
#endif //GLEW_PATH
#ifndef FREEGLUT_EXT_PATH
#   define FREEGLUT_EXT_PATH "GL/freeglut_ext.h"    // Optional (used only if glut.h comes from the freeglut library)
#endif //GLEW_PATH
#ifndef GLEW_PATH
#   define GLEW_PATH "GL/glew.h"    // Mandatory for Windows only
#endif //GLEW_PATH

#ifdef _WIN32
#	include "windows.h"
#	define USE_GLEW
#endif //_WIN32

#ifdef USE_GLEW
#	include GLEW_PATH
#else //USE_GLEW
#	define GL_GLEXT_PROTOTYPES
#endif //USE_GLEW

#include GLUT_PATH
#ifdef __FREEGLUT_STD_H__
#   include FREEGLUT_EXT_PATH
#endif //__FREEGLUT_STD_H__

// These derived definitions can't be touched [XSTR_MACRO(...) is used to insert an 'integer definition' between double quotes]
#define STR_MACRO(s) #s
#define XSTR_MACRO(s) STR_MACRO(s)
#define SHADOW_MAP_NUM_CASCADES_STRING XSTR_MACRO(SHADOW_MAP_NUM_CASCADES)

#include "helper_functions.h"   // please search this .c file for "Helper_":
                                // only very few of its functions are used.

#include <stdio.h>
#include <math.h>
#include <string.h>


// Config file handling: basically there's an .ini file next to the
// exe that you can tweak. (it's just an extra)
const char* ConfigFileName = PROGRAM_NAME".ini";
typedef struct {
    int fullscreen_width,fullscreen_height;
    int windowed_width,windowed_height;
    int fullscreen_enabled;
    int show_fps;
} Config;
void Config_Init(Config* c) {
    c->fullscreen_width=c->fullscreen_height=0;
    c->windowed_width=960;c->windowed_height=540;
    c->fullscreen_enabled=0;
    c->show_fps = 0;
}
int Config_Load(Config* c,const char* filePath)  {
    FILE* f = fopen(filePath, "rt");
    char ch='\0';char buf[256]="";
    size_t nread=0;
    int numParsedItem=0;
    if (!f)  return -1;
    while ((ch = fgetc(f)) !=EOF)    {
        buf[nread]=ch;
        nread++;
        if (nread>255) {
            nread=0;
            continue;
        }
        if (ch=='\n') {
           buf[nread]='\0';
           if (nread<2 || buf[0]=='[' || buf[0]=='#') {nread = 0;continue;}
           if (nread>2 && buf[0]=='/' && buf[1]=='/') {nread = 0;continue;}
           // Parse
           switch (numParsedItem)    {
               case 0:
               sscanf(buf, "%d %d", &c->fullscreen_width,&c->fullscreen_height);
               break;
               case 1:
               sscanf(buf, "%d %d", &c->windowed_width,&c->windowed_height);
               break;
               case 2:
               sscanf(buf, "%d", &c->fullscreen_enabled);
               break;
               case 4:
               sscanf(buf, "%d", &c->show_fps);
               break;
           }
           nread=0;
           ++numParsedItem;
        }
    }
    fclose(f);
    if (c->windowed_width<=0) c->windowed_width=720;
    if (c->windowed_height<=0) c->windowed_height=405;
    return 0;
}
int Config_Save(Config* c,const char* filePath)  {
    FILE* f = fopen(filePath, "wt");
    if (!f)  return -1;
    fprintf(f, "[Size In Fullscreen Mode (zero means desktop size)]\n%d %d\n",c->fullscreen_width,c->fullscreen_height);
    fprintf(f, "[Size In Windowed Mode]\n%d %d\n",c->windowed_width,c->windowed_height);
    fprintf(f, "[Fullscreen Enabled (0 or 1) (CTRL+RETURN)]\n%d\n", c->fullscreen_enabled);
    fprintf(f, "[Show FPS (0 or 1) (F2)]\n%d\n", c->show_fps);
    fprintf(f,"\n");
    fclose(f);
    return 0;
}

Config config;

// glut has a special fullscreen GameMode that you can toggle with CTRL+RETURN (not in WebGL)
int windowId = 0; 			// window Id when not in fullscreen mode
int gameModeWindowId = 0;	// window Id when in fullscreen mode

// Now we can start with our program

// camera data:
float targetPos[3];         // please set it in resetCamera()
float cameraYaw;            // please set it in resetCamera()
float cameraPitch;          // please set it in resetCamera()
float cameraDistance;       // please set it in resetCamera()
float cameraPos[3];         // Derived value (do not edit)
float vMatrix[16];          // view matrix
float cameraSpeed = 0.5f;   // When moving it

// light data
float lightYaw = M_PI*0.425f,lightPitch = M_PI*0.235f;    // must be copied to resetLight() too
float lightDirection[4] = {0,1,0,0};                    // Derived value (do not edit) [lightDirection[3]==0]

// pMatrix data:
float pMatrix[16];                  // projection matrix
const float pMatrixFovyDeg = 45.f;
const float pMatrixNearPlane = 0.5f;
const float pMatrixFarPlane = 20.0f;

// we calculate these in ResizeGL(...)
float gCascadeNearAndFarClippingPlanes[SHADOW_MAP_NUM_CASCADES+1];  // Array of the clipping planes of each cascade (gCascadeNearAndFarClippingPlanes[0]==pMatrixNearPlane and gCascadeNearAndFarClippingPlanes[SHADOW_MAP_NUM_CASCADES]==pMatrixFarPlane)
#ifdef USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
float gCascadePMatricesInv[16*SHADOW_MAP_NUM_CASCADES]; // One inverse pMatrix per cascade (calculated in ResizeGL(...))
#endif //USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE

float instantFrameTime = 16.2f;

// Optional (to speed up Helper_GlutDrawGeometry(...) a bit)
GLuint gDisplayListBase = 0;GLuint* pgDisplayListBase = &gDisplayListBase;  // Can be set to 0 as a fallback.


static const char* ShadowPassVertexShader[] = {
    "   void main() {\n"
    "       gl_Position = ftransform();\n"
    "   }\n"
};
static const char* ShadowPassFragmentShader[] = {
    "   void main() {\n"
    "       //gl_FragColor =  gl_Color;\n"
    "   }\n"
};
typedef struct {
    GLuint fbo;
    GLuint textureId;
    GLuint program;
} ShadowPass;

ShadowPass shadowPass;
void InitShadowPass(ShadowPass* sp)	{
    sp->program = Helper_LoadShaderProgramFromSource(*ShadowPassVertexShader,*ShadowPassFragmentShader);

    // create depth texture
    glGenTextures(1, &sp->textureId);
    glBindTexture(GL_TEXTURE_2D_ARRAY, sp->textureId);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, SHADOW_MAP_FILTER);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, SHADOW_MAP_FILTER);
#   ifndef __EMSCRIPTEN__
    //glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT32F, SHADOW_MAP_RESOLUTION, SHADOW_MAP_RESOLUTION, SHADOW_MAP_NUM_CASCADES, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    //glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT16, SHADOW_MAP_RESOLUTION, SHADOW_MAP_RESOLUTION, SHADOW_MAP_NUM_CASCADES, GL_FALSE, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0);
    glTexImage3D( GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_RESOLUTION, SHADOW_MAP_RESOLUTION, SHADOW_MAP_NUM_CASCADES, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, 0);
#   else //__EMSCRIPTEN__
    glTexImage3D( GL_TEXTURE_2D_ARRAY, 0, GL_DEPTH_COMPONENT, SHADOW_MAP_RESOLUTION, SHADOW_MAP_RESOLUTION, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_SHORT, 0);
#   undef SHADOW_MAP_CLAMP_MODE
#   define SHADOW_MAP_CLAMP_MODE GL_CLAMP_TO_EDGE
#   endif //__EMSCRIPTEN__
    if (SHADOW_MAP_CLAMP_MODE==GL_CLAMP_TO_BORDER)  {
        const GLfloat border[] = {1.0f,1.0f,1.0f,0.0f };
        glTexParameterfv(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_BORDER_COLOR, border);
    }
    glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, SHADOW_MAP_CLAMP_MODE );
    glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, SHADOW_MAP_CLAMP_MODE );
    glBindTexture(GL_TEXTURE_2D_ARRAY, 0);

    // create depth fbo
    glGenFramebuffers(1, &sp->fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, sp->fbo);
#   ifndef __EMSCRIPTEN__
    glDrawBuffer(GL_NONE); // Instruct openGL that we won't bind a color texture with the currently bound FBO
    glReadBuffer(GL_NONE);
#   endif //__EMSCRIPTEN__
    glFramebufferTexture(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,sp->textureId,0);

    {
        //Does the GPU support current FBO configuration?
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (status!=GL_FRAMEBUFFER_COMPLETE) printf("glCheckFramebufferStatus(...) FAILED for shadowPass.fbo. (status = %u)\n",status);
        // 36055 => FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT
    }
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}
void DestroyShadowPass(ShadowPass* sp)	{
    if (sp->program) {glDeleteProgram(sp->program);sp->program=0;}
    if (sp->fbo) {glDeleteBuffers(1,&sp->fbo);sp->fbo=0;}
    if (sp->textureId) {glDeleteTextures(1,&sp->textureId);}
}


static const char* DefaultPassVertexShader[] = {
    "varying vec4 v_diffuse;\n"
    "varying vec4 v_vertexModelViewSpace;\n"
    "varying float v_vertexModelViewSpaceDepth;\n"  // A bit redundant (we can calculate it in the fragment shader using v_vertexModelViewSpace)
    "\n"
    "void main()	{\n"
    "	gl_Position = ftransform();\n"
    "\n"
    "	vec3 normal = gl_NormalMatrix * gl_Normal;\n"
    "	vec3 lightVector = gl_LightSource[0].position.xyz\n;// - gl_Vertex.xyz;\n"
    "	float nxDir = max(0.0, dot(normal, lightVector));\n"
    "	v_diffuse = gl_LightSource[0].diffuse * nxDir; \n"
    "\n"
    "	gl_FrontColor = gl_Color;\n"
    "\n"
    "   v_vertexModelViewSpace = gl_ModelViewMatrix*gl_Vertex;\n"
    "   v_vertexModelViewSpaceDepth = -v_vertexModelViewSpace.z/v_vertexModelViewSpace.w;\n"  // Negative, because camera looks in the -z axis (and u_cascadeNearAndFarClippingPlanes[...] are all positive quantities)
    "}\n"
};
static const char* DefaultPassFragmentShader[] = {
#   ifdef VISUALIZE_CASCADE_SPLITS
    "#version 120\n"
    "#extension GL_EXT_texture_array : enable\n"
    "#extension GL_EXT_gpu_shader4 : enable\n"
    "#define NUM_CASCADES "SHADOW_MAP_NUM_CASCADES_STRING"\n"
    "const vec3 dbgSplitColors[4] = vec3[4](vec3(0.5,0.0,0.0),vec3(0.0,0.5,0.0),vec3(0.0,0.0,0.5),vec3(0.5,0.5,0.0));\n"
#   else //VISUALIZE_CASCADE_SPLITS
    "#extension GL_EXT_texture_array : enable\n"
    "#define NUM_CASCADES "SHADOW_MAP_NUM_CASCADES_STRING"\n"
#   endif //VISUALIZE_CASCADE_SPLITS
    "\n"
    "uniform sampler2DArray u_shadowMap;\n"
    "uniform vec2 u_shadowDarkening;\n" // .x = fDarkeningFactor [10.0-80.0], .y = min value clamp [0.0-1.0]
    "uniform float u_cascadeNearAndFarClippingPlanes[NUM_CASCADES+1];\n"
    "uniform mat4 u_biasedShadowMvpMatrix[NUM_CASCADES];\n" // Actually they are: (u_biasedShadowMvpMatrix[NUM_CASCADES] * vMatrixInverseCamera) please see the code.
    "\n"
    "varying vec4 v_diffuse;\n"
    "varying vec4 v_vertexModelViewSpace;\n"
    "varying float v_vertexModelViewSpaceDepth;\n"
    "\n"
    "void main() {\n"
    "	// Figure out which cascade to sample from\n"
    "   float cascadeIdxFloat = float(NUM_CASCADES-1);\n"
    "   for(int i=1;i<NUM_CASCADES;i++)  {\n"
    "       //if (v_vertexModelViewSpaceDepth < u_cascadeNearAndFarClippingPlanes[i]) cascadeIdxFloat-=1.0;\n"      // branch!
    "       cascadeIdxFloat-=max(sign(u_cascadeNearAndFarClippingPlanes[i] - v_vertexModelViewSpaceDepth), 0.0);\n" // branchless!
    "   }\n"
    "\n"
    "   int cascadeIdx = int(cascadeIdxFloat);\n"
    "   vec4 lightSpacePos = u_biasedShadowMvpMatrix[cascadeIdx]*v_vertexModelViewSpace;\n"   // There's a hidden vMatrixInverseCamera multiplication that removes the view component, moving the mMatrix from the camera space to the light space
    "   vec4 shadowCoordinateWdivide = lightSpacePos/lightSpacePos.w;\n"
    "   float shadowFactor = clamp(exp(u_shadowDarkening.x*(texture2DArray(u_shadowMap,vec3(shadowCoordinateWdivide.st,float(cascadeIdx))).r-shadowCoordinateWdivide.z)),u_shadowDarkening.y,1.0);\n"
    "\n"
#   ifdef VISUALIZE_CASCADE_SPLITS
    "   vec3 color = dbgSplitColors[cascadeIdx%4];\n"
#   else //VISUALIZE_CASCADE_SPLITS
    "   vec3 color = gl_Color.rgb;\n"
#   endif //VISUALIZE_CASCADE_SPLITS
    "\n"
    "	gl_FragColor = gl_LightSource[0].ambient + (v_diffuse * vec4(color*shadowFactor,1.0));\n"
    "}\n"
};


typedef struct {
    GLuint program;
    GLint uniform_location_biasedShadowMvpMatrix;
    GLint uniform_location_shadowMap;
    GLint uniform_location_shadowDarkening;
    GLint uniform_location_cascadeNearAndFarClippingPlanes;
} DefaultPass;

DefaultPass defaultPass;
void InitDefaultPass(DefaultPass* dp)	{
    dp->program = Helper_LoadShaderProgramFromSource(*DefaultPassVertexShader,*DefaultPassFragmentShader);
    dp->uniform_location_biasedShadowMvpMatrix = glGetUniformLocation(dp->program,"u_biasedShadowMvpMatrix");
    dp->uniform_location_shadowMap = glGetUniformLocation(dp->program,"u_shadowMap");
    dp->uniform_location_shadowDarkening = glGetUniformLocation(dp->program,"u_shadowDarkening");
    dp->uniform_location_cascadeNearAndFarClippingPlanes = glGetUniformLocation(dp->program,"u_cascadeNearAndFarClippingPlanes");

    glUseProgram(dp->program);
    glUniform1i(dp->uniform_location_shadowMap,0);
    glUniform2f(dp->uniform_location_shadowDarkening,80.0,0.45);	// Default values are (40.0f,0.75f) in [0-80] and [0-1]
    //glUniformMatrix4fv(dp->uniform_location_biasedShadowMvpMatrix, SHADOW_MAP_NUM_CASCADES /*only setting 1 matrix*/, GL_FALSE /*transpose?*/, Matrix);
    //glUniform1fv(dp->uniform_location_cascadeNearAndFarClippingPlanes,SHADOW_MAP_NUM_CASCADES,gCascadeNearAndFarClippingPlanes);
    glUseProgram(0);
}
void DestroyDefaultPass(DefaultPass* dp)	{
	if (dp->program) {glDeleteProgram(dp->program);dp->program=0;}
}


static const char* VisualizeTextureVertexShader[] = {
    "varying vec4 v_shadowCoord;\n"
    "\n"
    "   void main() {\n"
    "       gl_TexCoord[0] = gl_TextureMatrix[0] * gl_MultiTexCoord0;\n"
    "       gl_FrontColor = gl_Color;\n"
    "       gl_Position = ftransform();\n"
    "   }\n"
};
static const char* VisualizeTextureFragmentShader[] = {
    "//#define NUM_CASCADES "SHADOW_MAP_NUM_CASCADES_STRING"\n"
    "#extension GL_EXT_texture_array : enable\n"
    "\n"
    "uniform sampler2DArray u_shadowMap;\n"
    "uniform float u_textureArrayIndex;\n"
    "varying vec4 v_shadowCoord;\n"
    "\n"
    "   void main() {\n"
    "	vec4 shadowCoordinateWdivide = gl_TexCoord[0]/gl_TexCoord[0].w;\n"
    "       gl_FragColor = vec4(gl_Color.rbg*texture2DArray(u_shadowMap,vec3(shadowCoordinateWdivide.st,u_textureArrayIndex)).r,1.0);\n"
    "   }\n"
};
typedef struct {
    GLuint program;
    GLint uniform_location_shadowMap;
    GLint uniform_location_textureArrayIndex;
} VisualizeTexturePass;

VisualizeTexturePass visualizeTexturePass;
void InitVisualizeTexturePass(VisualizeTexturePass* vp)	{
    vp->program = Helper_LoadShaderProgramFromSource(*VisualizeTextureVertexShader,*VisualizeTextureFragmentShader);
    vp->uniform_location_shadowMap = glGetUniformLocation(vp->program,"u_shadowMap");
    vp->uniform_location_textureArrayIndex = glGetUniformLocation(vp->program,"u_textureArrayIndex");

    glUseProgram(vp->program);
    glUniform1i(vp->uniform_location_shadowMap,0);
    glUniform1f(vp->uniform_location_textureArrayIndex,0.f);
    glUseProgram(0);
}
void DestroyVisualizeTexturePass(VisualizeTexturePass* vp)	{
    if (vp->program) {glDeleteProgram(vp->program);vp->program=0;}
}



float current_width=0,current_height=0,current_aspect_ratio=1;  // Not sure when I've used these...
void ResizeGL(int w,int h) {
    current_width = (float) w;
    current_height = (float) h;
    if (current_height!=0) current_aspect_ratio = current_width/current_height;
    if (h>0)	{
        // We set our pMatrix
        Helper_Perspective(pMatrix,pMatrixFovyDeg,(float)w/(float)h,pMatrixNearPlane,pMatrixFarPlane);
        glMatrixMode(GL_PROJECTION);glLoadMatrixf(pMatrix);glMatrixMode(GL_MODELVIEW);

        // Here we calculate the splits (gCascadeNearAndFarClippingPlanes must contain SHADOW_MAP_NUM_CASCADES+1 elements) based on lambda, and camera near and far planes:
        Helper_GetCascadeNearAndFarClippingPlaneArray(gCascadeNearAndFarClippingPlanes,SHADOW_MAP_NUM_CASCADES,SHADOW_MAP_CASCADE_LAMBDA,pMatrixNearPlane,pMatrixFarPlane);
        // Here we update the uniform 'u_cascadeNearAndFarClippingPlanes' in the default pass shader program
        glUseProgram(defaultPass.program);
        glUniform1fv(defaultPass.uniform_location_cascadeNearAndFarClippingPlanes, SHADOW_MAP_NUM_CASCADES,gCascadeNearAndFarClippingPlanes);
        glUseProgram(0);

#       ifdef USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
        // Here we fill float gCascadePMatricesInv[16*SHADOW_MAP_NUM_CASCADES]: (basically the inverse of all the camera projection matrices, one per split)
        {
            int i;
            for (i=0;i<SHADOW_MAP_NUM_CASCADES;i++) {
                float* pMatrixInv = &gCascadePMatricesInv[16*i];
                Helper_Perspective(pMatrixInv,pMatrixFovyDeg,(float)w/(float)h,gCascadeNearAndFarClippingPlanes[i],gCascadeNearAndFarClippingPlanes[i+1]);
                Helper_InvertMatrix(pMatrixInv,pMatrixInv); // in-place operation (note tha we can't use Helper_InvertMatrixFast(...) here)
            }
        }
#       endif //USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
    }


    if (w>0 && h>0 && !config.fullscreen_enabled) {
        // On exiting we'll like to save these data back
        config.windowed_width=w;
        config.windowed_height=h;
    }

    glViewport(0,0,w,h);    // This is what people often call in ResizeGL()

}


void InitGL(void) {

    // These are important, but often overlooked OpenGL calls
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  // Otherwise transparent objects are not displayed correctly
    glClearColor(0.3f, 0.6f, 1.0f, 1.0f);
    glEnable(GL_TEXTURE_2D_ARRAY);    // Only needed for ffp, when VISUALIZE_DEPTH_TEXTURE is defined


	// ffp stuff
    glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);	
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);

	// New
    InitShadowPass(&shadowPass);
    InitDefaultPass(&defaultPass);
    InitVisualizeTexturePass(&visualizeTexturePass);

    // Please note that after InitGL(), this implementation calls ResizeGL(...,...).
    // If you copy/paste this code you can call it explicitly...
}

void DestroyGL() {
	// New
    DestroyShadowPass(&shadowPass);
    DestroyDefaultPass(&defaultPass);
    DestroyVisualizeTexturePass(&visualizeTexturePass);

    // 40 display lists are generated by Helper_GlutDrawGeometry(...) if pgDisplayListBase!=0
    if (pgDisplayListBase && *pgDisplayListBase) {glDeleteLists(*pgDisplayListBase,40);*pgDisplayListBase=0;}
}


void DrawGL(void) 
{	
    // We need to calculate the "instantFrameTime", because it's necessary to "dynamic_resolution.h"
    static unsigned begin = 0;
    static unsigned cur_time = 0;
    unsigned elapsed_time,delta_time;

    static float vMatrixInverse[16];
    static float lvpMatrices[SHADOW_MAP_NUM_CASCADES*16];                 // = light_pMatrix*light_vMatrix
    static float biasedShadowMvpMatrices[SHADOW_MAP_NUM_CASCADES*16];     // multiplied per vMatrixInverse

    float elapsedMs;float cosAlpha,sinAlpha;    // used to move objects around

    int i;

    if (begin==0) begin = glutGet(GLUT_ELAPSED_TIME);
    elapsed_time = glutGet(GLUT_ELAPSED_TIME) - begin;
    delta_time = elapsed_time - cur_time;
    instantFrameTime = (float)delta_time*0.001f;
    cur_time = elapsed_time;

    elapsedMs = (float)elapsed_time;
    cosAlpha = cos(elapsedMs*0.0005f);
    sinAlpha = sin(elapsedMs*0.00075f);

    // view Matrix
    Helper_LookAt(vMatrix,cameraPos[0],cameraPos[1],cameraPos[2],targetPos[0],targetPos[1],targetPos[2],0,1,0);    
	glLoadMatrixf(vMatrix);
    glLightfv(GL_LIGHT0,GL_POSITION,lightDirection);    // Important: the ffp must recalculate internally lightDirectionEyeSpace based on vMatrix [=> every frame]

    // view Matrix inverse (it's the camera matrix). Used twice below. So it's better to keep it here.
    Helper_InvertMatrixFast(vMatrixInverse,vMatrix);


    // Draw to Shadow Map------------------------------------------------------------------------------------------
    {
#       ifndef USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
        Helper_GetLightViewProjectionMatrices(lvpMatrices,gCascadeNearAndFarClippingPlanes,SHADOW_MAP_NUM_CASCADES,
                                              vMatrixInverse,
                                              pMatrixFovyDeg,current_aspect_ratio,
                                              lightDirection,1.0f/(float)SHADOW_MAP_RESOLUTION);
#       else //USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
        Helper_GetLightViewProjectionMatricesExtra(0,gCascadeNearAndFarClippingPlanes,SHADOW_MAP_NUM_CASCADES,
                                                   vMatrixInverse,
                                                   pMatrixFovyDeg,current_aspect_ratio,
                                                   lightDirection,1.0f/(float)SHADOW_MAP_RESOLUTION,
                                                   0,0,
                                                   gCascadePMatricesInv,    // Mandatory when we need to retrieve arguments that follow it
                                                   0,0,
                                                   lvpMatrices  // Technically this was provided as an 'lvpMatrices for optimal frustum culling usage' argument to be used in the 'Stable Shadow Mapping' case (but can be used to replace it too)
                                                   );
#       endif //USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE

        // Draw to shadow map texture
        glMatrixMode(GL_PROJECTION);glPushMatrix();glLoadIdentity();glMatrixMode(GL_MODELVIEW);        // We'll set the combined light view-projection matrix in GL_MODELVIEW (do you know that it's the same?)
        glBindFramebuffer(GL_FRAMEBUFFER, shadowPass.fbo);
        glViewport(0, 0, SHADOW_MAP_RESOLUTION,SHADOW_MAP_RESOLUTION);
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
        glCullFace(GL_FRONT);
        glEnable(GL_DEPTH_CLAMP);
        glUseProgram(shadowPass.program);
        for (i=0;i<SHADOW_MAP_NUM_CASCADES;i++) {
            glFramebufferTextureLayer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, shadowPass.fbo, 0, i);
            glClear(GL_DEPTH_BUFFER_BIT);   // Clears all the shadow map layer
            glPushMatrix();glLoadMatrixf(&lvpMatrices[i*16]); // we load both (light) projection and view matrices here (it's the same after all)
            Helper_GlutDrawGeometry(elapsedMs,cosAlpha,sinAlpha,targetPos,pgDisplayListBase);  // Done SHADOW_MAP_NUM_CASCADES times!
            glPopMatrix();
        }
        glUseProgram(0);
        glDisable(GL_DEPTH_CLAMP);
        glCullFace(GL_BACK);
        glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
        glBindFramebuffer(GL_FRAMEBUFFER,0);
        glMatrixMode(GL_PROJECTION);glPopMatrix();glMatrixMode(GL_MODELVIEW);

    }

    // Draw world
    {
        // biasedShadowMvpMatrix is used only in the DefaultPass:
        static float bias[16] = {0.5,0,0,0, 0,0.5,0,0,  0,0,0.5,0,    0.5,0.5,0.5,1}; // Moving from unit cube [-1,1] to [0,1]
        for (i=0;i<SHADOW_MAP_NUM_CASCADES;i++) {
            Helper_MultMatrix(&biasedShadowMvpMatrices[i*16],bias,&lvpMatrices[i*16]);
            Helper_MultMatrix(&biasedShadowMvpMatrices[i*16],&biasedShadowMvpMatrices[i*16],vMatrixInverse);  // We do this, so that when in the vs we multiply it with the camera mvMatrix, we get: biasedShadowMvpMatrix * mMatrix (using mMatrices directly in the shaders prevents the usage of double precision matrices: mvMatrices are good when converted to float to feed the shader, mMatrices are bad)
        }

        // Draw to world
        glViewport(0, 0, current_width,current_height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glBindTexture(GL_TEXTURE_2D_ARRAY,shadowPass.textureId);
        glUseProgram(defaultPass.program);
        glUniformMatrix4fv(defaultPass.uniform_location_biasedShadowMvpMatrix, SHADOW_MAP_NUM_CASCADES /*setting many matrices*/, GL_FALSE /*transpose?*/,biasedShadowMvpMatrices);
        Helper_GlutDrawGeometry(elapsedMs,cosAlpha,sinAlpha,targetPos,pgDisplayListBase);
        glUseProgram(0);
        glBindTexture(GL_TEXTURE_2D_ARRAY,0);
    }



    if (config.show_fps && instantFrameTime>0) {
        if ((elapsed_time/1000)%2==0)   {
            printf("FPS=%1.0f\n",1.f/instantFrameTime);fflush(stdout);
            config.show_fps=0;
        }
    }


#   ifdef VISUALIZE_DEPTH_TEXTURE
    {
        float startX=-1.f,endX=1.f;
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glDepthMask(GL_FALSE);

        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
	
		glColor3f(1,1,1);
		glDisable(GL_LIGHTING);
        glEnable(GL_BLEND);
        glBindTexture(GL_TEXTURE_2D_ARRAY,shadowPass.textureId);
        glColor4f(1,1,1,0.9f);
        glUseProgram(visualizeTexturePass.program);
        for (i=0;i<SHADOW_MAP_NUM_CASCADES;i++) {
            endX = startX+0.25f/current_aspect_ratio;

            glUniform1f(visualizeTexturePass.uniform_location_textureArrayIndex,(float)i);
            glBegin(GL_QUADS);
            glTexCoord2f(0,0);glVertex2f(startX,    0.75);
            glTexCoord2f(1,0);glVertex2f(endX, 0.75);
            glTexCoord2f(1,1);glVertex2f(endX, 1.0);
            glTexCoord2f(0,1);glVertex2f(startX,    1.0);
            glEnd();

            startX = endX;
        }

        glUseProgram(0);
        glBindTexture(GL_TEXTURE_2D_ARRAY,0);
        glDisable(GL_BLEND);
		glEnable(GL_LIGHTING);

        glPopMatrix();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_CULL_FACE);
        glDepthMask(GL_TRUE);
    }
#   endif //VISUALIZE_DEPTH_TEXTURE


}

static void GlutDestroyWindow(void);
static void GlutCreateWindow();

void GlutCloseWindow(void)  {Config_Save(&config,ConfigFileName);}

void GlutNormalKeys(unsigned char key, int x, int y) {
    const int mod = glutGetModifiers();
    switch (key) {
    case 27: 	// esc key
        Config_Save(&config,ConfigFileName);
        GlutDestroyWindow();
#		ifdef __FREEGLUT_STD_H__
        glutLeaveMainLoop();
#		else
        exit(0);
#		endif
        break;
    case 13:	// return key
    {
        if (mod&GLUT_ACTIVE_CTRL) {
            config.fullscreen_enabled = gameModeWindowId ? 0 : 1;
            GlutDestroyWindow();
            GlutCreateWindow();
        }
    }
        break;
    }

}

static void updateCameraPos() {
    const float distanceY = sin(cameraPitch)*cameraDistance;
    const float distanceXZ = cos(cameraPitch)*cameraDistance;
    cameraPos[0] = targetPos[0] + sin(cameraYaw)*distanceXZ;
    cameraPos[1] = targetPos[1] + distanceY;
    cameraPos[2] = targetPos[2] + cos(cameraYaw)*distanceXZ;
}

static void updateDirectionalLight() {
    const float distanceY = sin(lightPitch);
    const float distanceXZ = cos(lightPitch);
    lightDirection[0] = sin(lightYaw)*distanceXZ;
    lightDirection[1] = distanceY;
    lightDirection[2] = cos(lightYaw)*distanceXZ;
    Helper_Vector3Normalize(lightDirection);
    lightDirection[3]=0.f;
}

static void resetCamera() {
    // You can set the initial camera position here through:
    targetPos[0]=0; targetPos[1]=0; targetPos[2]=0; // The camera target point
    cameraYaw = 2*M_PI;                             // The camera rotation around the Y axis
    cameraPitch = M_PI*0.125f;                      // The camera rotation around the XZ plane
    cameraDistance = 5;                             // The distance between the camera position and the camera target point

    updateCameraPos();
}

static void resetLight() {
    lightYaw = M_PI*0.425f;
    lightPitch = M_PI*0.235f;
    updateDirectionalLight();
}

void GlutSpecialKeys(int key,int x,int y)
{
    const int mod = glutGetModifiers();
    if (!(mod&GLUT_ACTIVE_CTRL) && !(mod&GLUT_ACTIVE_SHIFT))	{
        switch (key) {
        case GLUT_KEY_LEFT:
        case GLUT_KEY_RIGHT:
            cameraYaw+= instantFrameTime*cameraSpeed*(key==GLUT_KEY_LEFT ? -4.0f : 4.0f);
            if (cameraYaw>M_PI) cameraYaw-=2*M_PI;
            else if (cameraYaw<=-M_PI) cameraYaw+=2*M_PI;
            updateCameraPos();		break;
        case GLUT_KEY_UP:
        case GLUT_KEY_DOWN:
            cameraPitch+= instantFrameTime*cameraSpeed*(key==GLUT_KEY_UP ? 2.f : -2.f);
            if (cameraPitch>M_PI-0.001f) cameraPitch=M_PI-0.001f;
            else if (cameraPitch<-M_PI*0.05f) cameraPitch=-M_PI*0.05f;
            updateCameraPos();
            break;
        case GLUT_KEY_PAGE_UP:
        case GLUT_KEY_PAGE_DOWN:
            cameraDistance+= instantFrameTime*cameraSpeed*(key==GLUT_KEY_PAGE_DOWN ? 25.0f : -25.0f);
            if (cameraDistance<1.f) cameraDistance=1.f;
            updateCameraPos();
            break;
        case GLUT_KEY_F1:
        case GLUT_KEY_F2:
            config.show_fps = !config.show_fps;
            //printf("showFPS: %s.\n",config.show_fps?"ON":"OFF");fflush(stdout);
            break;
        case GLUT_KEY_HOME:
            // Reset the camera
            resetCamera();
            break;
        }
    }
    else if (mod&GLUT_ACTIVE_CTRL) {
        switch (key) {
        case GLUT_KEY_LEFT:
        case GLUT_KEY_RIGHT:
        case GLUT_KEY_UP:
        case GLUT_KEY_DOWN:
        {
            // Here we move targetPos and cameraPos at the same time

            // We must find a pivot relative to the camera here (ignoring Y)
            float forward[3] = {targetPos[0]-cameraPos[0],0,targetPos[2]-cameraPos[2]};
            float up[3] = {0,1,0};
            float left[3];

            Helper_Vector3Normalize(forward);
            Helper_Vector3Cross(left,up,forward);
            {
                float delta[3] = {0,0,0};int i;
                if (key==GLUT_KEY_LEFT || key==GLUT_KEY_RIGHT) {
                    float amount = instantFrameTime*cameraSpeed*(key==GLUT_KEY_RIGHT ? -25.0f : 25.0f);
                    for (i=0;i<3;i++) delta[i]+=amount*left[i];
                }
                else {
                    float amount = instantFrameTime*cameraSpeed*(key==GLUT_KEY_DOWN ? -25.0f : 25.0f);
                    for ( i=0;i<3;i++) delta[i]+=amount*forward[i];
                }
                for ( i=0;i<3;i++) {
                    targetPos[i]+=delta[i];
                    cameraPos[i]+=delta[i];
                }
            }
        }
        break;
        case GLUT_KEY_PAGE_UP:
        case GLUT_KEY_PAGE_DOWN:
            // We use world space coords here.
            targetPos[1]+= instantFrameTime*cameraSpeed*(key==GLUT_KEY_PAGE_DOWN ? -25.0f : 25.0f);
            if (targetPos[1]<-50.f) targetPos[1]=-50.f;
            else if (targetPos[1]>500.f) targetPos[1]=500.f;
            updateCameraPos();
        break;
        }
    }
    else if (mod&GLUT_ACTIVE_SHIFT)	{
        switch (key) {
        case GLUT_KEY_LEFT:
        case GLUT_KEY_RIGHT:
            lightYaw+= instantFrameTime*cameraSpeed*(key==GLUT_KEY_LEFT ? -4.0f : 4.0f);
            if (lightYaw>M_PI) lightYaw-=2*M_PI;
            else if (lightYaw<=-M_PI) lightYaw+=2*M_PI;
            updateDirectionalLight();
            break;
        case GLUT_KEY_UP:
        case GLUT_KEY_DOWN:
        case GLUT_KEY_PAGE_UP:
        case GLUT_KEY_PAGE_DOWN:
            lightPitch+= instantFrameTime*cameraSpeed*( (key==GLUT_KEY_UP || key==GLUT_KEY_PAGE_UP) ? 2.f : -2.f);
            if (lightPitch>M_PI-0.001f) lightPitch=M_PI-0.001f;
            else if (lightPitch<-M_PI*0.05f) lightPitch=-M_PI*0.05f;
            updateDirectionalLight();
            break;
        case GLUT_KEY_HOME:
            // Reset the light
            resetLight();
            break;
        }
    }
}

void GlutMouse(int a,int b,int c,int d) {

}

// Note that we have used GlutFakeDrawGL() so that at startup
// the calling order is: InitGL(),ResizeGL(...),DrawGL()
// Also note that glutSwapBuffers() must NOT be called inside DrawGL()
static void GlutDrawGL(void)		{DrawGL();glutSwapBuffers();}
static void GlutIdle(void)			{glutPostRedisplay();}
static void GlutFakeDrawGL(void) 	{glutDisplayFunc(GlutDrawGL);}
void GlutDestroyWindow(void) {
    if (gameModeWindowId || windowId)	{
        DestroyGL();

        if (gameModeWindowId) {
            glutLeaveGameMode();
            gameModeWindowId = 0;
        }
        if (windowId) {
            glutDestroyWindow(windowId);
            windowId=0;
        }
    }
}
void GlutCreateWindow() {
    GlutDestroyWindow();
    if (config.fullscreen_enabled)	{
        char gms[16]="";
        if (config.fullscreen_width>0 && config.fullscreen_height>0)	{
            sprintf(gms,"%dx%d:32",config.fullscreen_width,config.fullscreen_height);
            glutGameModeString(gms);
            if (glutGameModeGet (GLUT_GAME_MODE_POSSIBLE)) gameModeWindowId = glutEnterGameMode();
            else config.fullscreen_width=config.fullscreen_height=0;
        }
        if (gameModeWindowId==0)	{
            const int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
            const int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
            sprintf(gms,"%dx%d:32",screenWidth,screenHeight);
            glutGameModeString(gms);
            if (glutGameModeGet (GLUT_GAME_MODE_POSSIBLE)) gameModeWindowId = glutEnterGameMode();
        }
    }
    if (!gameModeWindowId) {
        char windowTitle[1024] = PROGRAM_NAME".c\t";
#       ifdef  USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
        strcat(windowTitle,"[Unstable]\t");
#       endif  //USE_UNSTABLE_SHADOW_MAPPING_TECHNIQUE
        strcat(windowTitle,"("XSTR_MACRO(SHADOW_MAP_RESOLUTION)"x"SHADOW_MAP_NUM_CASCADES_STRING")");
        config.fullscreen_enabled = 0;
        glutInitWindowPosition(100,100);
        glutInitWindowSize(config.windowed_width,config.windowed_height);
        windowId = glutCreateWindow(windowTitle);
    }

    glutKeyboardFunc(GlutNormalKeys);
    glutSpecialFunc(GlutSpecialKeys);
    glutMouseFunc(GlutMouse);
    glutIdleFunc(GlutIdle);
    glutReshapeFunc(ResizeGL);
    glutDisplayFunc(GlutFakeDrawGL);
#   ifdef __FREEGLUT_STD_H__
    glutWMCloseFunc(GlutCloseWindow);
#   endif //__FREEGLUT_STD_H__

#ifdef USE_GLEW
    {
        GLenum err = glewInit();
        if( GLEW_OK != err ) {
            fprintf(stderr, "Error initializing GLEW: %s\n", glewGetErrorString(err) );
            return;
        }
    }
#endif //USE_GLEW

    InitGL();

}


int main(int argc, char** argv)
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    //glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
#ifdef __FREEGLUT_STD_H__
    glutSetOption ( GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION ) ;
#endif //__FREEGLUT_STD_H__

    Config_Init(&config);
    Config_Load(&config,ConfigFileName);

    GlutCreateWindow();

    //OpenGL info
    printf("\nGL Vendor: %s\n", glGetString( GL_VENDOR ));
    printf("GL Renderer : %s\n", glGetString( GL_RENDERER ));
    printf("GL Version (string) : %s\n",  glGetString( GL_VERSION ));
    printf("GLSL Version : %s\n", glGetString( GL_SHADING_LANGUAGE_VERSION ));
    //printf("GL Extensions:\n%s\n",(char *) glGetString(GL_EXTENSIONS));
    {int maxTextureSize=0;glGetIntegerv(GL_MAX_TEXTURE_SIZE, &maxTextureSize);printf("Max Texture Size: %d\n", maxTextureSize);}
    // Many GPUs return 16384 as their max size.
    // For 4 bytes per pixel: 16384*16384*4 bytes = 1 Gb
    // For 4 float per pixel: 16384*16384*4 floats (a floating point texture) = 4Gb.
    // Many GPUs don't even have 4Gb of memory

    printf("\nKEYS:\n");
    printf("AROW KEYS + PAGE_UP/PAGE_DOWN:\tmove camera (optionally with CTRL down)\n");
    printf("HOME KEY:\t\t\treset camera\n");
    printf("ARROW KEYS + SHIFT:\tmove directional light\n");
    printf("CTRL+RETURN:\t\ttoggle fullscreen on/off\n");
    printf("F2:\t\t\tdisplay FPS\n");
    printf("\n");

    resetCamera();  // Mandatory
    resetLight();   // Mandatory

    glutMainLoop();


    return 0;
}










