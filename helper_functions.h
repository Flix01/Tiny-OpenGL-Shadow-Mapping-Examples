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

#ifndef HELPER_FUNCTIONS_H_
#define HELPER_FUNCTIONS_H_
#include <math.h>
#include <stdio.h>	// printf

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HELPER_FUNCTIONS_VERSION
#   define HELPER_FUNCTIONS_VERSION 1.0
#endif //VERSION


/* The __restrict and __restrict__ keywords are recognized in both C, at all language levels, and C++, at LANGLVL(EXTENDED).*/
//#ifdef NO_RESTRICT  // please define it globally if the keyword __restrict is not present
#   ifndef __restrict   
#       define __restrict /*no-op*/
#   endif
//#endif //NO_RESTRICT


#ifndef MATRIX_USE_DOUBLE_PRECISION
typedef float  hloat;       // short form of hELPER_FUNCTIONS_FLoat
#else
typedef double hloat;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_HALF_PI
#define M_HALF_PI (M_PI/2.0)
#endif
#ifndef M_PIOVER180
#define M_PIOVER180 (3.14159265358979323846/180.0)
#endif
#ifndef M_180OVERPI
#define M_180OVERPI (180.0/3.14159265358979323846)
#endif


static __inline hloat Helper_Round(hloat number)	{return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);}
static __inline void Helper_IdentityMatrix(hloat* __restrict result16) {
    hloat* m = result16;
    m[0]=m[5]=m[10]=m[15]=1;
    m[1]=m[2]=m[3]=m[4]=m[6]=m[7]=m[8]=m[9]=m[11]=m[12]=m[13]=m[14]=0;
}
static __inline void Helper_CopyMatrix(hloat* __restrict dst16,const hloat* __restrict src16) {
    int i;for (i=0;i<16;i++) dst16[i]=src16[i];
}
static __inline void Helper_MultMatrix(hloat* __restrict result16,const hloat* __restrict ml16,const hloat* __restrict mr16) {
    int i,j,j4;
    if (result16==ml16) {
        hloat ML16[16];Helper_CopyMatrix(ML16,ml16);
        Helper_MultMatrix(result16,ML16,mr16);
        return;
    }
    else if (result16==mr16) {
        hloat MR16[16];Helper_CopyMatrix(MR16,mr16);
        Helper_MultMatrix(result16,ml16,MR16);
        return;
    }
    for(i = 0; i < 4; i++) {
        for(j = 0; j < 4; j++) {
            j4 = 4*j;
            result16[i+j4] =
                ml16[i]    * mr16[0+j4] +
                ml16[i+4]  * mr16[1+j4] +
                ml16[i+8]  * mr16[2+j4] +
                ml16[i+12] * mr16[3+j4];
        }
    }
}
static __inline void Helper_MatrixMulDir(const hloat* __restrict m16,hloat* __restrict dirOut3,const hloat dirX,hloat dirY,hloat dirZ) {
    //hloat w;
    dirOut3[0] = dirX*m16[0] + dirY*m16[4] + dirZ*m16[8];
    dirOut3[1] = dirX*m16[1] + dirY*m16[5] + dirZ*m16[9];
    dirOut3[2] = dirX*m16[2] + dirY*m16[6] + dirZ*m16[10];
    //w          = dirX*m16[3] + dirY*m16[7] + dirZ*m16[11]; // + m[15] ?
    //if (w!=0 && w!=1) {dirOut3[0]/=w;dirOut3[1]/=w;dirOut3[2]/=w;}
}
static __inline void Helper_MatrixMulPos(const hloat* __restrict m16,hloat* __restrict posOut3,const hloat posX,hloat posY,hloat posZ) {
    hloat w;
    posOut3[0] = posX*m16[0] + posY*m16[4] + posZ*m16[8] + m16[12];
    posOut3[1] = posX*m16[1] + posY*m16[5] + posZ*m16[9] + m16[13];
    posOut3[2] = posX*m16[2] + posY*m16[6] + posZ*m16[10]+ m16[14];
    w          = posX*m16[3] + posY*m16[7] + posZ*m16[11]+ m16[15];
    if (w!=0 && w!=1) {posOut3[0]/=w;posOut3[1]/=w;posOut3[2]/=w;}
}
static __inline void Helper_NormalizePlane(hloat* __restrict p4) {
    hloat len = p4[0]*p4[0]+p4[1]*p4[1]+p4[2]*p4[2];
    if (len!=0) {
        int i;len = sqrt(len);
        for (i=0;i<4;i++) p4[i]/=len;
    }
}
static __inline hloat Helper_Vector3Dot(const hloat* a3,const hloat* b3) {return a3[0]*b3[0]+a3[1]*b3[1]+a3[2]*b3[2];}
static __inline void Helper_Vector3Normalize(hloat* __restrict v3) {
    hloat len = Helper_Vector3Dot(v3,v3);int i;
    if (len!=0) {len = sqrt(len);for (i=0;i<3;i++) v3[i]/=len;}
}
static __inline void Helper_Vector3Cross(hloat* __restrict vOut3,const hloat* __restrict a3,const hloat* __restrict b3) {
    vOut3[0] =	a3[1] * b3[2] - a3[2] * b3[1];
    vOut3[1] =	a3[2] * b3[0] - a3[0] * b3[2];
    vOut3[2] =	a3[0] * b3[1] - a3[1] * b3[0];
}
static __inline hloat Helper_Vector3DistSquared(const hloat* __restrict a3,const hloat* __restrict b3) {
    const hloat rv[3] = {b3[0]-a3[0],b3[1]-a3[1],b3[2]-a3[2]};
    return rv[0]*rv[0] + rv[1]*rv[1] + rv[2]*rv[2];
}
static __inline hloat Helper_Vector3Dist(const hloat* __restrict a3,const hloat* __restrict b3) {
    const hloat res = Helper_Vector3DistSquared(a3,b3);
    return res!=0 ? sqrt(res) : 0;
}


static __inline void Helper_ConvertMatrixd2f16(float* __restrict result16,const double* __restrict m16) {int i;for(i = 0; i < 16; i++) result16[i]=(float)m16[i];}
static __inline void Helper_ConvertMatrixf2d16(double* __restrict result16,const float* __restrict m16) {int i;for(i = 0; i < 16; i++) result16[i]=(double)m16[i];}
static __inline void Helper_ConvertMatrixd2f9(float* __restrict result9,const double* __restrict m9) {int i;for(i = 0; i < 9; i++) result9[i]=(float)m9[i];}
static __inline void Helper_ConvertMatrixf2d9(double* __restrict result9,const float* __restrict m9) {int i;for(i = 0; i < 9; i++) result9[i]=(double)m9[i];}


__inline static void Helper_GlUniformMatrix4v(GLint location,GLsizei count,GLboolean transpose,const hloat* value) {
    const float* fvalue = NULL;
#   ifndef MATRIX_USE_DOUBLE_PRECISION
    fvalue = value;
#   else
    float val[16];Helper_ConvertMatrixd2f16(val,value);fvalue=val;
#   endif
    glUniformMatrix4fv(location,count,transpose,fvalue);
}
__inline static void Helper_GlUniformMatrix3v(GLint location,GLsizei count,GLboolean transpose,const hloat* value) {
    const float* fvalue = NULL;
#   ifndef MATRIX_USE_DOUBLE_PRECISION
    fvalue = value;
#   else
    float val[9];Helper_ConvertMatrixd2f9(val,value);fvalue=val;
#   endif
    glUniformMatrix3fv(location,count,transpose,fvalue);
}
__inline static void Helper_GlUniform3v(GLint location,GLsizei count,const hloat* value) {
    const float* fvalue = NULL;
#   ifndef MATRIX_USE_DOUBLE_PRECISION
    fvalue = value;
#   else
    const float val[3] = {(float)value[0],(float)value[1],(float)value[2]};fvalue = val;
#   endif
    glUniform3fv(location,count,fvalue);
}


static __inline void Helper_LookAt(hloat* __restrict mOut16,hloat eyeX,hloat eyeY,hloat eyeZ,hloat centerX,hloat centerY,hloat centerZ,hloat upX,hloat upY,hloat upZ)    {
    hloat* m = mOut16;
    const hloat eps = 0.0001;

    hloat F[3] = {eyeX-centerX,eyeY-centerY,eyeZ-centerZ};
    hloat length = F[0]*F[0]+F[1]*F[1]+F[2]*F[2];	// length2 now
    hloat up[3] = {upX,upY,upZ};

    hloat S[3] = {up[1]*F[2]-up[2]*F[1],up[2]*F[0]-up[0]*F[2],up[0]*F[1]-up[1]*F[0]};
    hloat U[3] = {F[1]*S[2]-F[2]*S[1],F[2]*S[0]-F[0]*S[2],F[0]*S[1]-F[1]*S[0]};

    if (length==0) length = eps;
    length = sqrt(length);
    F[0]/=length;F[1]/=length;F[2]/=length;

    length = S[0]*S[0]+S[1]*S[1]+S[2]*S[2];if (length==0) length = eps;
    length = sqrt(length);
    S[0]/=length;S[1]/=length;S[2]/=length;

    length = U[0]*U[0]+U[1]*U[1]+U[2]*U[2];if (length==0) length = eps;
    length = sqrt(length);
    U[0]/=length;U[1]/=length;U[2]/=length;

    /*
                S0	S1	S2  0		1	0	0	-ex
                U0	U1	U2	0   *   0	1	0	-ey
                F0	F1	F2  0		0	0	1	-ez
                0	0	0	1		0	0	0	1

            */
    m[0] = S[0];
    m[1] = U[0];
    m[2] = F[0];
    m[3]= 0;

    m[4] = S[1];
    m[5] = U[1];
    m[6] = F[1];
    m[7]= 0;

    m[8] = S[2];
    m[9] = U[2];
    m[10]= F[2];
    m[11]= 0;

    m[12] = -S[0]*eyeX -S[1]*eyeY -S[2]*eyeZ;
    m[13] = -U[0]*eyeX -U[1]*eyeY -U[2]*eyeZ;
    m[14]= -F[0]*eyeX -F[1]*eyeY -F[2]*eyeZ;
    m[15]= 1;
}
static __inline void Helper_Perspective(hloat* __restrict mOut16,hloat degfovy,hloat aspect, hloat zNear, hloat zFar) {
    hloat* res = mOut16;
    const hloat eps = 0.0001;
    hloat f = 1.f/tan(degfovy*3.14159265358979323846/360.0); //cotg(degfovy/2)
    hloat Dfn = (zFar-zNear);
    if (Dfn==0) {zFar+=eps;zNear-=eps;Dfn=zFar-zNear;}
    if (aspect==0) aspect = 1.f;

    res[0]  = f/aspect;
    res[1]  = 0;
    res[2]  = 0;
    res[3]  = 0;

    res[4]  = 0;
    res[5]  = f;
    res[6]  = 0;
    res[7] = 0;

    res[8]  = 0;
    res[9]  = 0;
    res[10] = -(zFar+zNear)/Dfn;
    res[11] = -1;

    res[12]  = 0;
    res[13]  = 0;
    res[14] = -2.f*zFar*zNear/Dfn;
    res[15] = 0;
}
static __inline void Helper_Ortho(hloat* __restrict mOut16,hloat left,hloat right, hloat bottom, hloat top,hloat nearVal,hloat farVal) {
    hloat* res = mOut16;
    const hloat eps = 0.0001;
    hloat Drl = (right-left);
    hloat Dtb = (top-bottom);
    hloat Dfn = (farVal-nearVal);
    if (Drl==0) {right+=eps;left-=eps;Drl=right-left;}
    if (Dtb==0) {top+=eps;bottom-=eps;Dtb=top-bottom;}
    if (Dfn==0) {farVal+=eps;nearVal-=eps;Dfn=farVal-nearVal;}

    res[0]  = 2.f/Drl;
    res[1]  = 0;
    res[2]  = 0;
    res[3] = 0;

    res[4]  = 0;
    res[5]  = 2.f/Dtb;
    res[6]  = 0;
    res[7] = 0;

    res[8]  = 0;
    res[9]  = 0;
    res[10] = -2.f/Dfn;
    res[11] = 0;

    res[12]  = -(right+left)/Drl;
    res[13]  = -(top+bottom)/Dtb;
    res[14] = (farVal+nearVal)/Dfn;
    res[15] = 1;
}
static __inline int Helper_InvertMatrix(hloat* __restrict mOut16,const hloat* __restrict m16)	{
    const hloat* m = m16;
    hloat* n = mOut16;

    hloat m00 = m[0],  m10 = m[1],  m20 = m[2],  m30 = m[3];
    hloat m01 = m[4],  m11 = m[5],  m21 = m[6],  m31 = m[7];
    hloat m02 = m[8],  m12 = m[9],  m22 = m[10], m32 = m[11];
    hloat m03 = m[12], m13 = m[13], m23 = m[14], m33 = m[15];

    hloat v0 = m20 * m31 - m21 * m30;
    hloat v1 = m20 * m32 - m22 * m30;
    hloat v2 = m20 * m33 - m23 * m30;
    hloat v3 = m21 * m32 - m22 * m31;
    hloat v4 = m21 * m33 - m23 * m31;
    hloat v5 = m22 * m33 - m23 * m32;

    hloat t00 = + (v5 * m11 - v4 * m12 + v3 * m13);
    hloat t10 = - (v5 * m10 - v2 * m12 + v1 * m13);
    hloat t20 = + (v4 * m10 - v2 * m11 + v0 * m13);
    hloat t30 = - (v3 * m10 - v1 * m11 + v0 * m12);

    hloat det = (t00 * m00 + t10 * m01 + t20 * m02 + t30 * m03);
    if (det==0) return 0;
    {
        hloat invDet = 1 / det;

        hloat d00 = t00 * invDet;
        hloat d10 = t10 * invDet;
        hloat d20 = t20 * invDet;
        hloat d30 = t30 * invDet;

        hloat d01 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
        hloat d11 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
        hloat d21 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
        hloat d31 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

        v0 = m10 * m31 - m11 * m30;
        v1 = m10 * m32 - m12 * m30;
        v2 = m10 * m33 - m13 * m30;
        v3 = m11 * m32 - m12 * m31;
        v4 = m11 * m33 - m13 * m31;
        v5 = m12 * m33 - m13 * m32;
        {
            hloat d02 = + (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
            hloat d12 = - (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
            hloat d22 = + (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
            hloat d32 = - (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

            v0 = m21 * m10 - m20 * m11;
            v1 = m22 * m10 - m20 * m12;
            v2 = m23 * m10 - m20 * m13;
            v3 = m22 * m11 - m21 * m12;
            v4 = m23 * m11 - m21 * m13;
            v5 = m23 * m12 - m22 * m13;
            {
                hloat d03 = - (v5 * m01 - v4 * m02 + v3 * m03) * invDet;
                hloat d13 = + (v5 * m00 - v2 * m02 + v1 * m03) * invDet;
                hloat d23 = - (v4 * m00 - v2 * m01 + v0 * m03) * invDet;
                hloat d33 = + (v3 * m00 - v1 * m01 + v0 * m02) * invDet;

                n[0] =d00; n[1] =d10; n[2] =d20; n[3] =d30;
                n[4] =d01; n[5] =d11; n[6] =d21; n[7] =d31;
                n[8] =d02; n[9] =d12; n[10]=d22; n[11]=d32;
                n[12]=d03; n[13]=d13; n[14]=d23; n[15]=d33;
            }
        }
    }
    return 1;
}
static __inline void Helper_InvertMatrixFast(hloat* __restrict mOut16,const hloat* __restrict m16)	{
    // It works only for translation + rotation, and only
    // when rotation can be represented by an unit quaternion
    // scaling is discarded
    const hloat* m = m16;
    hloat* n = mOut16;
    const hloat T[3] = {-m[12],-m[13],-m[14]};
    hloat w;
    // Step 1. Transpose the 3x3 submatrix
    n[3]=m[3];n[7]=m[7];n[11]=m[11];n[15]=m[15];
    n[0]=m[0];n[1]=m[4];n[2]=m[8];
    n[4]=m[1];n[5]=m[5];n[6]=m[9];
    n[8]=m[2];n[9]=m[6];n[10]=m[10];
    // Step2. Adjust translation
    n[12]=T[0]*n[0] + T[1]*n[4] +T[2]*n[8];
    n[13]=T[0]*n[1] + T[1]*n[5] +T[2]*n[9];
    n[14]=T[0]*n[2] + T[1]*n[6] +T[2]*n[10];
    w    =T[0]*n[3] + T[1]*n[7] +T[2]*n[11];
    if (w!=0 && w!=1) {n[12]/=w;n[13]/=w;n[14]/=w;} // These last 2 lines are not strictly necessary AFAIK
}
static __inline void Helper_GetFrustumPlaneEquations(hloat planeEquationsOut[6][4],const hloat* __restrict vpMatrix16,int normalizePlanes)   {
    // ax+by+cz+d=0 [xl,xr,yb,yt,zn,zf],normalizePlanes=0 -> no normalization
    hloat m00 = vpMatrix16[0], m01 = vpMatrix16[4], m02 = vpMatrix16[8],  m03 = vpMatrix16[12];
    hloat m10 = vpMatrix16[1], m11 = vpMatrix16[5], m12 = vpMatrix16[9],  m13 = vpMatrix16[13];
    hloat m20 = vpMatrix16[2], m21 = vpMatrix16[6], m22 = vpMatrix16[10], m23 = vpMatrix16[14];
    hloat m30 = vpMatrix16[3], m31 = vpMatrix16[7], m32 = vpMatrix16[11], m33 = vpMatrix16[15];
    hloat* p = NULL;
    p = &planeEquationsOut[0][0];   // Left
    p[0] = m30+m00; p[1] = m31+m01; p[2] = m32+m02; p[3] = m33+m03;if (normalizePlanes) Helper_NormalizePlane(p);
    p = &planeEquationsOut[1][0];   // Right
    p[0] = m30-m00; p[1] = m31-m01; p[2] = m32-m02; p[3] = m33-m03;if (normalizePlanes) Helper_NormalizePlane(p);
    p = &planeEquationsOut[2][0];   // Bottom
    p[0] = m30+m10; p[1] = m31+m11; p[2] = m32+m12; p[3] = m33+m13;if (normalizePlanes) Helper_NormalizePlane(p);
    p = &planeEquationsOut[3][0];   // Top
    p[0] = m30-m10; p[1] = m31-m11; p[2] = m32-m12; p[3] = m33-m13;if (normalizePlanes) Helper_NormalizePlane(p);
    p = &planeEquationsOut[4][0];   // Near
    p[0] = m30+m20; p[1] = m31+m21; p[2] = m32+m22; p[3] = m33+m23;if (normalizePlanes) Helper_NormalizePlane(p);
    p = &planeEquationsOut[5][0];   // Far
    p[0] = m30-m20; p[1] = m31-m21; p[2] = m32-m22; p[3] = m33-m23;if (normalizePlanes) Helper_NormalizePlane(p);
}
static __inline void Helper_GetFrustumPoints(hloat frustumPoints[8][4],const hloat* __restrict vpMatrixInverse16)    {
    const hloat v[8][4] = {{-1, -1, -1, 1},{-1,  1, -1, 1},{ 1,  1, -1, 1},{ 1, -1, -1, 1},{-1, -1, 1, 1},{-1,  1, 1, 1},{ 1,  1, 1, 1},{ 1, -1, 1, 1}};
    int i;for (i = 0; i < 8; i++) {
        Helper_MatrixMulPos(vpMatrixInverse16,frustumPoints[i],v[i][0],v[i][1],v[i][2]);
        frustumPoints[i][3]=1;
    }
}
static __inline void Helper_GetFrustumAabbCenterAndHalfExtents(hloat* __restrict frustumCenterOut3,hloat* __restrict frustumHalfExtentsOut3,const hloat frustumPoints[8][4])    {
    hloat vmin[3] = {frustumPoints[0][0],frustumPoints[0][1],frustumPoints[0][2]};
    hloat vmax[3] = {vmin[0],vmin[1],vmin[2]};
    int i,j;
    for (i = 1; i < 8; i++) {
        for (j = 0; j < 3; j++) {
            if      (vmin[j] > frustumPoints[i][j]) vmin[j] = frustumPoints[i][j];
            else if (vmax[j] < frustumPoints[i][j]) vmax[j] = frustumPoints[i][j];
        }
    }
    for (j = 0; j < 3; j++) {
        if (frustumCenterOut3)      frustumCenterOut3[j] = (vmin[j]+vmax[j])*0.5;
        if (frustumHalfExtentsOut3) frustumHalfExtentsOut3[j] = (vmax[j]-vmin[j])*0.5;
    }
}

static __inline void Helper_GetLightViewProjectionMatrix(hloat* __restrict lvpMatrixOut16,
                                                          const hloat*  __restrict cameraPosition3,
                                                          const hloat*  __restrict cameraForwardDirection3,
                                                          hloat cameraNearClippingPlane,hloat cameraFarClippingPlane,hloat cameraFovyDeg,hloat cameraAspectRatio,
                                                          const hloat*  __restrict normalizedLightDirection3, hloat texelIncrement)  {
    hloat frustumCenter[3] = {0,0,0};hloat radius = 0;
    hloat lpMatrix[16],lvMatrix[16];
    int i;

    // Get frustumCenter and radius
    hloat frustumCenterDistance;
    hloat tanFovDiagonalSquared = tan(cameraFovyDeg*3.14159265358979323846/360.0); // 0.5*M_PI/180.0
    const hloat halfNearFarClippingPlane = 0.5*(cameraFarClippingPlane+cameraNearClippingPlane);
    tanFovDiagonalSquared*=tanFovDiagonalSquared;
    tanFovDiagonalSquared*=(1.0+cameraAspectRatio*cameraAspectRatio);
    frustumCenterDistance = halfNearFarClippingPlane*(1.0+tanFovDiagonalSquared);
    if (frustumCenterDistance > cameraFarClippingPlane) frustumCenterDistance = cameraFarClippingPlane;
    for (i=0;i<3;i++) frustumCenter[i] = cameraPosition3[i]+cameraForwardDirection3[i]*frustumCenterDistance;
    radius = sqrt( (tanFovDiagonalSquared*cameraFarClippingPlane*cameraFarClippingPlane) + (cameraFarClippingPlane-frustumCenterDistance)*(cameraFarClippingPlane-frustumCenterDistance));
    //fprintf(stderr,"radius=%1.4f frustumCenterDistance=%1.4f nearPlane=%1.4f farPlane = %1.4f\n",radius,frustumCenterDistance,cameraNearClippingPlane,cameraFarClippingPlane);

    // For people trying to save texture space it's:  halfNearFarClippingPlane <= frustumCenterDistance <= cameraFarClippingPlane
    // When frustumCenterDistance == cameraFarClippingPlane, then frustumCenter is on the far clip plane (and half the texture space gets wasted).
    // when frustumCenterDistance == halfNearFarClippingPlane, then we're using an ortho projection matrix, and frustumCenter is in the middle of the near and far plane (no texture space gets wasted).
    // in all the other cases the space wasted can go from zero to half texture

    // Shadow swimming happens when: 1) camera translates; 2) camera rotates; 3) objects move or rotate
    // AFAIK Shadow swimming (3) can't be fixed in any way
    if (texelIncrement>0)   radius = ceil(radius/texelIncrement)*texelIncrement;      // This 'should' fix Shadow swimming (1)  [Not sure code is correct!]

    // Get light matrices
    Helper_Ortho(lpMatrix,-radius,radius,-radius,radius,0,-2.0*radius);
    Helper_LookAt(lvMatrix,
            frustumCenter[0]-normalizedLightDirection3[0]*radius,   // eye[0]
            frustumCenter[1]-normalizedLightDirection3[1]*radius,   // eye[1]
            frustumCenter[2]-normalizedLightDirection3[2]*radius,   // eye[2]
            frustumCenter[0],frustumCenter[1],frustumCenter[2],     // target
            0,1,0                                                   // up (people that cares about wasted texture space can probably change it)
            );
    // Get output
    Helper_MultMatrix(lvpMatrixOut16,lpMatrix,lvMatrix);

    // This 'should' fix Shadow swimming (2) [Not sure code is correct!]
    if (texelIncrement>0)   {
        hloat shadowOrigin[4]   = {0,0,0,1};
        hloat roundedOrigin[4]  = {0,0,0,0};
        hloat roundOffset[4]    = {0,0,0,0};
        hloat texelCoefficient = texelIncrement*2.0;
        Helper_MatrixMulPos(lvpMatrixOut16,shadowOrigin,shadowOrigin[0],shadowOrigin[1],shadowOrigin[2]);
        for (i = 0; i < 2; i++) {// Or i<3 ?
            shadowOrigin[i]/= texelCoefficient;
            roundedOrigin[i] = Helper_Round(shadowOrigin[i]);
            roundOffset[i] = roundedOrigin[i] - shadowOrigin[i];
            roundOffset[i]*=  texelCoefficient;
        }

        lvpMatrixOut16[12]+= roundOffset[0];
        lvpMatrixOut16[13]+= roundOffset[1];
    }

    // Debug stuff
    //fprintf(stderr,"radius=%1.5f frustumCenter={%1.5f,%1.5f,%1.5f}\n",radius,frustumCenter[0],frustumCenter[1],frustumCenter[2]);
}



static __inline  void Helper_GetCascadeFarPlaneArray(hloat*  __restrict cascadeFarPlanesArrayOut,int numCascades,hloat lambda,hloat cameraNearClippingPlane,hloat cameraFarClippingPlane) {
    int i;
    for (i=0;i<numCascades;i++) {
        hloat p = (float)(i+1)/(float)(numCascades);
        hloat logValue = cameraNearClippingPlane * pow(cameraFarClippingPlane/cameraNearClippingPlane, p);
        hloat uniformValue = cameraNearClippingPlane + (cameraFarClippingPlane-cameraNearClippingPlane)*p;
        hloat d = lambda*(logValue-uniformValue)+uniformValue;
        // Far clip planes in (0,1]
        hloat cascadeSplit = (d-cameraNearClippingPlane)/(cameraFarClippingPlane-cameraNearClippingPlane);
        //if (cascadeFarPlanesInClipSpaceArrayOut) cascadeFarPlanesInClipSpaceArrayOut[i] = cascadeSplit;
        // Store cascadeFarPlanesArrayOut in OpenGL units in (cameraNearClippingPlane,cameraFarClippingPlane]
        cascadeFarPlanesArrayOut[i] = (cameraNearClippingPlane + cascadeSplit * (cameraFarClippingPlane - cameraNearClippingPlane));
        //fprintf(stderr,"cascadeSplits[%d] = %1.4f\tcascadeFarPlanesArrayOut[%d] = %1.4f\n",i,cascadeSplit,i,cascadeFarPlanesArrayOut[i]);
    }
}
static __inline  void Helper_GetLightViewProjectionMatrices(hloat* __restrict lvpMatricesOut16,const hloat*  __restrict cascadeFarPlanesArray,int numCascades,
                                                          const hloat*  __restrict cameraPosition3,
                                                          const hloat*  __restrict cameraForwardDirection3,
                                                          hloat cameraNearClippingPlane,hloat cameraFarClippingPlane,hloat cameraFovyDeg,hloat cameraAspectRatio,
                                                          const hloat*  __restrict normalizedLightDirection3, hloat texelIncrement
                                                          //,hloat* __restrict optionalCascadeSphereCenterPlanesOut,hloat* __restrict optionalCascadeSphereRadiiSquaredOut
                                                          //, const hloat* __restrict optionalCameraVMatrix16
                                                          )  {
    int cascadeIterator,i;
    // fovx = 2 * atan(aspect*tan(fovy/2))
    // tan(fovx/2) = aspect*tan(fovy/2)
    // tan(fovd/2)*tan(fovd/2) = (aspect*aspect+1.0)*tan(fovy/2)*tan(fovy/2))
    hloat tanFovDiagonalSquared = tan(cameraFovyDeg*3.14159265358979323846/360.0);  // 0.5*M_PI/180.0
    //hloat maxRadius = 0.f;
    tanFovDiagonalSquared*=tanFovDiagonalSquared;
    tanFovDiagonalSquared*=(1.0+cameraAspectRatio*cameraAspectRatio);
    for (cascadeIterator = numCascades-1; cascadeIterator >= 0; --cascadeIterator) {
        hloat frustumNearClippingPlane = cascadeIterator == 0 ? cameraNearClippingPlane : cascadeFarPlanesArray[cascadeIterator - 1];
        hloat frustumFarClippingPlane = cascadeFarPlanesArray[cascadeIterator];
        hloat frustumCenter[3] = {0,0,0};
        hloat radius = 0.0;//,radius2 = 0.0;
        hloat lpMatrix[16],lvMatrix[16];
        hloat* lvpMatrixOut16 = &lvpMatricesOut16[16*cascadeIterator];

        // Get frustumCenter and radius
        hloat frustumCenterDistance;
        const hloat halfNearFarClippingPlane = 0.5*(frustumFarClippingPlane+frustumNearClippingPlane);

        // Mine ------------------------------------------------------------------------------------------------
        frustumCenterDistance = halfNearFarClippingPlane*(1.0+tanFovDiagonalSquared);
        if (frustumCenterDistance > frustumFarClippingPlane) {frustumCenterDistance = frustumFarClippingPlane;} // (***)
        radius = (tanFovDiagonalSquared*frustumFarClippingPlane*frustumFarClippingPlane) + (frustumFarClippingPlane-frustumCenterDistance)*(frustumFarClippingPlane-frustumCenterDistance);
        //radius2 = (tanFovDiagonalSquared*frustumNearClippingPlane*frustumNearClippingPlane) + (frustumNearClippingPlane-frustumCenterDistance)*(frustumNearClippingPlane-frustumCenterDistance);if (radius<radius2) radius = radius2;

        //if (optionalCascadeSphereCenterPlanesOut) optionalCascadeSphereCenterPlanesOut[cascadeIterator] = frustumCenterDistance;
        //if (optionalCascadeSphereRadiiSquaredOut) optionalCascadeSphereRadiiSquaredOut[cascadeIterator] = radius;

        radius = sqrt(radius);  // Mandatory

        // ------------------------------------------------------------------------------------------------------
        /*
        // Or: https://lxjk.github.io/2017/04/15/Calculate-Minimal-Bounding-Sphere-of-Frustum.html (it uses fovx)-
        if (tanFovDiagonalSquared>=(frustumFarClippingPlane-frustumNearClippingPlane)/(frustumFarClippingPlane+frustumNearClippingPlane)) {frustumCenterDistance=frustumFarClippingPlane;radius=frustumFarClippingPlane*sqrt(tanFovDiagonalSquared);}
        else {
            frustumCenterDistance = halfNearFarClippingPlane*(1.0+tanFovDiagonalSquared);
            radius = 0.5 * sqrt((frustumFarClippingPlane-frustumNearClippingPlane)*(frustumFarClippingPlane-frustumNearClippingPlane)+
                                2.0*(frustumFarClippingPlane*frustumFarClippingPlane+frustumNearClippingPlane*frustumNearClippingPlane)*tanFovDiagonalSquared+
                                (frustumFarClippingPlane+frustumNearClippingPlane)*(frustumFarClippingPlane+frustumNearClippingPlane)*tanFovDiagonalSquared*tanFovDiagonalSquared
                                );
        }
        //--------------------------------------------------------------------------------------------------------
        */

        for (i=0;i<3;i++) frustumCenter[i] = cameraPosition3[i]+cameraForwardDirection3[i]*frustumCenterDistance;
        //fprintf(stderr,"cascade[%d] radius=%1.4f frustumCenterDistance=%1.4f nearPlane=%1.4f farPlane = %1.4f fovy = %1.0f cameraAspectRatio = %1.3f\n",cascadeIterator,radius,frustumCenterDistance,frustumNearClippingPlane,frustumFarClippingPlane,cameraFovyDeg,cameraAspectRatio);

        /*if (optionalCameraVMatrix16) {
            // Dbg only [To check that all frustum corners distances from frustumCenter should be equal if we comment out line (***) above]
            int j;
            hloat distance = 0,radius_cmp = 0;hloat distances[3] = {0,0,0};
            hloat cameraVPMatrixInv[16],camerapMatrix[16];hloat fustumPoints[8][4];
            Helper_Perspective(camerapMatrix,cameraFovyDeg,cameraAspectRatio,frustumNearClippingPlane,frustumFarClippingPlane);
            Helper_MultMatrix(cameraVPMatrixInv,camerapMatrix,optionalCameraVMatrix16);
            Helper_InvertMatrix(cameraVPMatrixInv,cameraVPMatrixInv);
            Helper_GetFrustumPoints(fustumPoints,cameraVPMatrixInv);

            // Get Max distance (for sure in perspective cameras, it's when i=4 or 5 or 6 or 7: the same max distance should be found)
            fprintf(stderr,"cascade[%d] corner distances:\t",cascadeIterator);
            for (i = 0; i < 8; i++) {
                for (j=0;j<3;j++) {
                    distances[j] = (fustumPoints[i][j] - frustumCenter[j]);
                    distances[j]*=distances[j];
                }
                distance = distances[0]+distances[1]+distances[2];  // squared
                if (radius_cmp < distance) radius_cmp = distance;           // squared
                fprintf(stderr,"%1.3f ",sqrt(distance));
            }
            radius_cmp = sqrt(radius_cmp);
            fprintf(stderr,"\tradius_cmp=%1.6f (radius=%1.6f)\n",radius_cmp,radius);
        }*/

        // Shadow swimming happens when: 1) camera translates; 2) camera rotates; 3) objects move or rotate
        // AFAIK Shadow swimming (3) can't be fixed in any way
        if (texelIncrement>0)   radius = ceil(radius/texelIncrement)*texelIncrement;      // This 'should' fix Shadow swimming (1)  [Not sure code is correct!]
        //if (cascadeIterator == numCascades-1) maxRadius = radius;

        // Get light matrices
        Helper_Ortho(lpMatrix,-radius,radius,-radius,radius,cascadeIterator==0 ? (2.0*radius) : 0.0,-2.0*radius); // maybe here we can use farVal (last arg) as radius (or maxRadius).
                                                                                                                   // also note that nearVal should be 0, but when radius is small (first cascades) it can produce artifacts (don't ask me why).
        Helper_LookAt(lvMatrix,
                frustumCenter[0]-normalizedLightDirection3[0]*radius,    // eye[0]
                frustumCenter[1]-normalizedLightDirection3[1]*radius,    // eye[1]
                frustumCenter[2]-normalizedLightDirection3[2]*radius,    // eye[2]
                frustumCenter[0],frustumCenter[1],frustumCenter[2],      // target
                0,1,0                                                    // up
                );
        // Get output
        Helper_MultMatrix(lvpMatrixOut16,lpMatrix,lvMatrix);

        // This 'should' fix Shadow swimming (2) [Not sure code is correct!]
        if (texelIncrement>0)   {
            hloat shadowOrigin[4]   = {0,0,0,1};
            hloat roundedOrigin[4]  = {0,0,0,0};
            hloat roundOffset[4]    = {0,0,0,0};
            hloat texelCoefficient = texelIncrement*2.0;
            Helper_MatrixMulPos(lvpMatrixOut16,shadowOrigin,shadowOrigin[0],shadowOrigin[1],shadowOrigin[2]);    // Or MultDir ?
            for (i = 0; i < 2; i++) {// Or i<3 ?
                shadowOrigin[i]/= texelCoefficient;
                roundedOrigin[i] = Helper_Round(shadowOrigin[i]);
                roundOffset[i] = roundedOrigin[i] - shadowOrigin[i];
                roundOffset[i]*=  texelCoefficient;
            }

            lvpMatrixOut16[12]+= roundOffset[0];
            lvpMatrixOut16[13]+= roundOffset[1];
        }
    }
}
static __inline void Helper_Min3(hloat* __restrict res3,const hloat* a3,const hloat* b3) {
    int i;for (i=0;i<3;i++) res3[i]=a3[i]<b3[i]?a3[i]:b3[i];
}
static __inline void Helper_Max3(hloat* __restrict res3,const hloat* a3,const hloat* b3) {
    int i;for (i=0;i<3;i++) res3[i]=a3[i]>b3[i]?a3[i]:b3[i];
}
// It "should" performs AABB test. mfMatrix16 is the matrix M so that: F*M = mvpMatrix (F being the matrix used to extract the frustum planes). We can use: F=pMatrix and M=mvMatrix, or: F=vpMatrix and M=mMatrix.
static __inline int Helper_IsVisible(const hloat frustumPlanes[6][4],const hloat*__restrict mfMatrix16,hloat aabbMinX,hloat aabbMinY,hloat aabbMinZ,hloat aabbMaxX,hloat aabbMaxY,hloat aabbMaxZ) {
    // It "should" performs AABB test. mfMatrix16 is the matrix M so that: F*M = mvpMatrix (F being the matrix used to extract the frustum planes). Here we use: F=pMatrix and M=mvMatrix, but it could be: F=vpMatrix and M=mMatrix too.
    int i;
    // AABB transformation by: http://dev.theomader.com/transform-bounding-boxes/
    const hloat* m = mfMatrix16;
    const hloat a[9] = {m[0]*aabbMinX,m[1]*aabbMinX,    m[2]*aabbMinX, m[4]*aabbMinY,m[5]*aabbMinY,m[6]*aabbMinY,   m[8]*aabbMinZ,m[9]*aabbMinZ,m[10]*aabbMinZ};
    const hloat b[9] = {m[0]*aabbMaxX,m[1]*aabbMaxX,    m[2]*aabbMaxX, m[4]*aabbMaxY,m[5]*aabbMaxY,m[6]*aabbMaxY,   m[8]*aabbMaxZ,m[9]*aabbMaxZ,m[10]*aabbMaxZ};
    hloat buf[18];
    Helper_Min3(&buf[0], &a[0],&b[0]);
    Helper_Min3(&buf[3], &a[3],&b[3]);
    Helper_Min3(&buf[6], &a[6],&b[6]);
    Helper_Max3(&buf[9], &a[0],&b[0]);
    Helper_Max3(&buf[12],&a[3],&b[3]);
    Helper_Max3(&buf[15],&a[6],&b[6]);
    {
        const hloat aabb[6] = {
            buf[0]+buf[ 3]+buf[ 6]+m[12], buf[ 1]+buf[ 4]+buf[ 7]+m[13], buf[ 2]+buf[ 5]+buf[ 8]+m[14],
            buf[9]+buf[12]+buf[15]+m[12], buf[10]+buf[13]+buf[16]+m[13], buf[11]+buf[14]+buf[17]+m[14]
        };

        // End AABB transformation

        // FAST VERSION
        {
            for(i=0; i < 6; i++) {
                const hloat *pl = &frustumPlanes[i][0];
                const int p[3] = {3*(int)(pl[0]>0.0f),3*(int)(pl[1]>0.0f),3*(int)(pl[2]>0.0f)};   // p[j] = 0 or 3
                const hloat dp = pl[0]*aabb[p[0]] + pl[1]*aabb[p[1]+1] + pl[2]*aabb[p[2]+2] + pl[3];
                if (dp < 0) return 0;
            }
        }

        // MUCH SLOWER VERSION
        /*{
        const hloat aabb0=aabb[0],aabb1=aabb[1],aabb2=aabb[2];
        const hloat aabb3=aabb[3],aabb4=aabb[4],aabb5=aabb[5];
        for(i=0; i < 6; i++)    {
            int sum = 0;
            const hloat *pl = &frustumPlanes[i][0];
            const hloat plx=pl[0],ply=pl[1],plz=pl[2],plw=pl[3];
            sum += (plx*aabb0+ply*aabb1+plz*aabb2+plw)<0?1:0;
            sum += (plx*aabb3+ply*aabb1+plz*aabb2+plw)<0?1:0;
            sum += (plx*aabb0+ply*aabb4+plz*aabb2+plw)<0?1:0;
            sum += (plx*aabb3+ply*aabb4+plz*aabb2+plw)<0?1:0;
            sum += (plx*aabb0+ply*aabb1+plz*aabb5+plw)<0?1:0;
            sum += (plx*aabb3+ply*aabb1+plz*aabb5+plw)<0?1:0;
            sum += (plx*aabb0+ply*aabb4+plz*aabb5+plw)<0?1:0;
            sum += (plx*aabb3+ply*aabb4+plz*aabb5+plw)<0?1:0;
            if (sum==8) return 0;
        }
    }*/

        // Furthermore we still have a lot of false positives in both cases
    }
    return 1;
}



#ifndef NO_HELPER_FUNCTIONS_OPENGL
// Loading shader function
static __inline GLhandleARB Helper_LoadShader(const char* buffer, const unsigned int type)
{
    GLhandleARB handle;
    const GLcharARB* files[1];

    // shader Compilation variable
    GLint result;				// Compilation code result
    GLint errorLoglength ;
    char* errorLogText;
    GLsizei actualErrorLogLength;

    handle = glCreateShader(type);
    if (!handle)
    {
        //We have failed creating the vertex shader object.
        printf("Failed creating vertex shader object.\n");
        exit(0);
    }

    files[0] = (const GLcharARB*)buffer;
    glShaderSource(
                handle, //The handle to our shader
                1, //The number of files.
                files, //An array of const char * data, which represents the source code of theshaders
                NULL);

    glCompileShader(handle);

    //Compilation checking.
    glGetShaderiv(handle, GL_COMPILE_STATUS, &result);

    // If an error was detected.
    if (!result)
    {
        //We failed to compile.
        printf("Shader failed compilation.\n");

        //Attempt to get the length of our error log.
        glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &errorLoglength);

        //Create a buffer to read compilation error message
        errorLogText =(char*) malloc(sizeof(char) * errorLoglength);

        //Used to get the final length of the log.
        glGetShaderInfoLog(handle, errorLoglength, &actualErrorLogLength, errorLogText);

        // Display errors.
        printf("%s\n",errorLogText);

        // Free the buffer malloced earlier
        free(errorLogText);
    }

    return handle;
}
static __inline GLuint Helper_LoadShaderProgramFromSourceWithGeometryShader(const char* vs,const char* gs,const char* fs)	{
    // shader Compilation variable
    GLint result;				// Compilation code result
    GLint errorLoglength ;
    char* errorLogText;
    GLsizei actualErrorLogLength;

    GLhandleARB vertexShaderHandle = 0;
    GLhandleARB geometryShaderHandle = 0;
    GLhandleARB fragmentShaderHandle = 0;
    GLuint programId = 0;

    if (vs) vertexShaderHandle   = Helper_LoadShader(vs,GL_VERTEX_SHADER);
    if (gs) geometryShaderHandle   = Helper_LoadShader(vs,GL_GEOMETRY_SHADER);
    if (fs) fragmentShaderHandle = Helper_LoadShader(fs,GL_FRAGMENT_SHADER);
    if (!vertexShaderHandle || !fragmentShaderHandle) return 0;

    programId = glCreateProgram();

    glAttachShader(programId,vertexShaderHandle);
    if (geometryShaderHandle) glAttachShader(programId,geometryShaderHandle);
    glAttachShader(programId,fragmentShaderHandle);
    glLinkProgram(programId);

    //Link checking.
    glGetProgramiv(programId, GL_LINK_STATUS, &result);

    // If an error was detected.
    if (!result)
    {
        //We failed to compile.
        printf("Program failed to link.\n");

        //Attempt to get the length of our error log.
        glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &errorLoglength);

        //Create a buffer to read compilation error message
        errorLogText =(char*) malloc(sizeof(char) * errorLoglength);

        //Used to get the final length of the log.
        glGetProgramInfoLog(programId, errorLoglength, &actualErrorLogLength, errorLogText);

        // Display errors.
        printf("%s\n",errorLogText);

        // Free the buffer malloced earlier
        free(errorLogText);
    }

    glDeleteShader(vertexShaderHandle);
    if (geometryShaderHandle) glDeleteShader(geometryShaderHandle);
    glDeleteShader(fragmentShaderHandle);

    return programId;
}
static __inline GLuint Helper_LoadShaderProgramFromSource(const char* vs,const char* fs)	{return Helper_LoadShaderProgramFromSourceWithGeometryShader(vs,0,fs);}


#if (defined(__glut_h__) || defined(__GLUT_H__) || defined(__FREEGLUT_STD_H__) || defined(__OPENGLUT_STD_H__))
// elapsedMs: from the start of the application
// cosAlpha,sinAlpha: some changing values used to move two objects
// cameraTargetPosition3 is an optional array of 3 floats containing the camera target position (to display a 3D pivot)
// pOptionalDisplayListBase: when not 0, 40 display lists are generated, created and used
//                          please free them at cleanup using glDeleteLists(*pOptionalDisplayListBase,40);
void Helper_GlutDrawGeometry(float elapsedMs, float cosAlpha, float sinAlpha,const float* cameraTargetPosition3,GLuint* pOptionalDisplayListBase) {
    int i=0,j=0,curDisplayListIndex=0;
    const int mustGenerateLists = (pOptionalDisplayListBase && *pOptionalDisplayListBase==0) ? 1 : 0;
    const int mustDrawGeometry = (!pOptionalDisplayListBase || *pOptionalDisplayListBase==0) ? 1 : 0;
    const int mustCallLists = (pOptionalDisplayListBase && *pOptionalDisplayListBase!=0) ? 1 : 0;

    /*static int lastMustDrawGeometry = 0;  // DBG
    if (lastMustDrawGeometry!=mustDrawGeometry) {
        lastMustDrawGeometry=mustDrawGeometry;
        fprintf(stderr,"mustDrawGeometry=%d mustGenerateLists=%d mustCallLists=%d\n",mustDrawGeometry,mustGenerateLists,mustCallLists);
    }*/

    if (mustGenerateLists) *pOptionalDisplayListBase = glGenLists(40);

    if (pOptionalDisplayListBase)   {
        if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
        else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
        ++curDisplayListIndex;
    }
    if (mustDrawGeometry)   {

        // ground
        glColor3f(0.2,0.4,0.2);
        glPushMatrix();

        glTranslatef(0,-0.5,0);

        glPushMatrix();
        glScalef(11.f,0.5f,14.f);
        glutSolidCube(1.0);
        glPopMatrix();

        glPopMatrix();

        // columns
        glColor3f(0.35,0.35,0.35);
        for (j=0;j<2;j++)	{
            for (i=0;i<=10;i++)	{
                glPushMatrix();

                glTranslatef(4.75-j*2.0,2.7+0.31,-5.f+(float)i*1.0);
                glPushMatrix();
                glRotatef(90,1,0,0);

                glPushMatrix();
                glScalef(0.2f,0.2f,2.8f);

                // central part
                glutSolidCylinder(0.5,1.0,8,8);

                // higher part
                glutSolidCone(0.8f,0.1f,8,8);
                glTranslatef(0.f,0.f,-0.025f);
                glutSolidCylinder(0.8,0.025,8,8);

                // lower part
                glTranslatef(0.f,0.f,1.05);
                glFrontFace(GL_CW);
                glutSolidCone(0.8f,-0.1f,8,8);
                glFrontFace(GL_CCW);
                glTranslatef(0.f,0.f,0.0f);
                glutSolidCylinder(0.8,0.025,8,8);

                glPopMatrix();


                glPopMatrix();

                glPopMatrix();
            }
        }


        // column plane under roof
        glColor3f(0.8,0.8,0.8);
        glPushMatrix();
        glTranslatef(3.75,3.16,0.f);
        glScalef(2.5f,0.155f,10.75f);
        glutSolidCube(1.0);
        glPopMatrix();

        // column roof
        glColor3f(0.4,0.0,0.0);
        glPushMatrix();
        glTranslatef(3.75,3.02+0.48,-10.75*0.5);
        glScalef(0.825f,0.3f,1.f);
        glRotatef(90,0,0,1);
        glutSolidCylinder(1.75,10.75,3,3);
        glPopMatrix();

        // column base
        glColor3f(0.2,0.2,0.2);
        glPushMatrix();

        glTranslatef(3.75,-0.01f,0.f);
        glScalef(2.8f,0.155f,10.5f);
        glutSolidCube(1.0);

        glTranslatef(0,-1,0);
        glScalef(1.15f,1,1.05f);
        glutSolidCube(1.0);

        glPopMatrix();

        // cube
        glColor3f(0.5,0.5,1.0);
        glPushMatrix();
        glTranslatef(-2,0.3,0.25);
        glScalef(0.5f,1.5f,0.5f);
        glutSolidCube(0.75);
        glPopMatrix();

        //#   define DBG_BIG_OCCLUDING_WALL
#       ifdef DBG_BIG_OCCLUDING_WALL
        // big occluding wall
        glColor3f(0.5,0.1,0.1);
        glPushMatrix();

        glTranslatef(10,0,0);

        glPushMatrix();
        glScalef(0.25f,2.f*pMatrixFarPlane,2.f*pMatrixFarPlane);
        glutSolidCube(1.0);
        glPopMatrix();

        glPopMatrix();
#       endif // DBG_BIG_OCCLUDING_WALL

    }
    if (mustGenerateLists) glEndList();

    // sphere
    glColor3f(0.8,0.8,0);
    glPushMatrix();

    glTranslatef(-1+2.5*cosAlpha,0.25,-1+sinAlpha);

    glPushMatrix();
    glRotatef(-elapsedMs*0.05f,0,1,0);

    if (pOptionalDisplayListBase)   {
        if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
        else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
        ++curDisplayListIndex;
    }
    if (mustDrawGeometry) glutSolidSphere(0.5,16,16);
    if (mustGenerateLists) glEndList();

    glPopMatrix();

    glPopMatrix();

    // tours
    glColor3f(0.4,0.4,0.8);
    glPushMatrix();

    glTranslatef(-0.5+2.5*cosAlpha,0.5,2-sinAlpha);

    glPushMatrix();
    glRotatef(elapsedMs*0.05f,0,1,0);
    glScalef(1.f,1.f,1.f);

    if (pOptionalDisplayListBase)   {
        if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
        else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
        ++curDisplayListIndex;
    }
    if (mustDrawGeometry) glutSolidTorus(0.25,0.5,16,16);
    if (mustGenerateLists) glEndList();

    glPopMatrix();

    glPopMatrix();

    // teapot
    glColor3f(0.4,0.2,0.0);
    glPushMatrix();

    glTranslatef(-0.4,0.1,-4);

    glPushMatrix();
    glRotatef(elapsedMs*0.1f,0,1,0);
    glScalef(1.f,1.f,1.f);
    if (pOptionalDisplayListBase)   {
        if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
        else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
        ++curDisplayListIndex;
    }
    if (mustDrawGeometry) {glFrontFace(GL_CW);glutSolidTeapot(0.5);glFrontFace(GL_CCW);}
    if (mustGenerateLists) glEndList();

    glPopMatrix();

    glPopMatrix();

    if (cameraTargetPosition3)  {
        // camera target pos:
        // Center
        glColor3f(0,0,0);
        glPushMatrix();
        glTranslatef(cameraTargetPosition3[0],cameraTargetPosition3[1],cameraTargetPosition3[2]);
        glPushMatrix();
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidSphere(0.04,8,8);
        if (mustGenerateLists) glEndList();

        // X Axis
        glPushMatrix();
        glColor3f(1,0,0);
        glRotatef(90,0,1,0);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCylinder(0.04,0.25,8,8);
        if (mustGenerateLists) glEndList();
        glTranslatef(0,0,0.25);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCone(0.06,0.1,8,8);
        if (mustGenerateLists) glEndList();
        glPopMatrix();

        // Y Axis
        glPushMatrix();
        glColor3f(0,1,0);
        glRotatef(-90,1,0,0);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCylinder(0.04,0.25,8,8);
        if (mustGenerateLists) glEndList();
        glTranslatef(0,0,0.25);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCone(0.06,0.1,8,8);
        if (mustGenerateLists) glEndList();
        glPopMatrix();

        // Z Axis
        glPushMatrix();
        glColor3f(0,0,1);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCylinder(0.04,0.25,8,8);
        if (mustGenerateLists) glEndList();
        glTranslatef(0,0,0.25);
        if (pOptionalDisplayListBase)   {
            if (mustCallLists) glCallList(*pOptionalDisplayListBase+curDisplayListIndex);
            else if (mustGenerateLists) glNewList(*pOptionalDisplayListBase+curDisplayListIndex,GL_COMPILE_AND_EXECUTE);
            ++curDisplayListIndex;
        }
        if (mustDrawGeometry) glutSolidCone(0.06,0.1,8,8);
        if (mustGenerateLists) glEndList();
        glPopMatrix();

        glPopMatrix();
        glPopMatrix();
        // End camera target position
    }
}
#endif // glut_h

#endif //NO_HELPER_FUNCTIONS_OPENGL

#ifdef __cplusplus
}   // extern "C"
#endif

#endif //HELPER_FUNCTIONS_H_

