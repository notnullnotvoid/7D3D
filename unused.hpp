

inline Mat4 transpose(Mat4 m) {
    return {
        m.m00, m.m10, m.m20, m.m30,
        m.m01, m.m11, m.m21, m.m31,
        m.m02, m.m12, m.m22, m.m32,
        m.m03, m.m13, m.m23, m.m33,
    };
}

//this is ripped straight from the MESA library, and adapted to work with our codebase
//it doesn't check the determinant, so it will return +-inf values if determinant is 0
inline Mat4 invert(Mat4 m) {
    Mat4 inv;
    inv.m00 =  m.m11 * m.m22 * m.m33 -
               m.m11 * m.m23 * m.m32 -
               m.m21 * m.m12 * m.m33 +
               m.m21 * m.m13 * m.m32 +
               m.m31 * m.m12 * m.m23 -
               m.m31 * m.m13 * m.m22;
    inv.m10 = -m.m10 * m.m22 * m.m33 +
               m.m10 * m.m23 * m.m32 +
               m.m20 * m.m12 * m.m33 -
               m.m20 * m.m13 * m.m32 -
               m.m30 * m.m12 * m.m23 +
               m.m30 * m.m13 * m.m22;
    inv.m20 =  m.m10 * m.m21 * m.m33 -
               m.m10 * m.m23 * m.m31 -
               m.m20 * m.m11 * m.m33 +
               m.m20 * m.m13 * m.m31 +
               m.m30 * m.m11 * m.m23 -
               m.m30 * m.m13 * m.m21;
    inv.m30 = -m.m10 * m.m21 * m.m32 +
               m.m10 * m.m22 * m.m31 +
               m.m20 * m.m11 * m.m32 -
               m.m20 * m.m12 * m.m31 -
               m.m30 * m.m11 * m.m22 +
               m.m30 * m.m12 * m.m21;
    inv.m01 = -m.m01 * m.m22 * m.m33 +
               m.m01 * m.m23 * m.m32 +
               m.m21 * m.m02 * m.m33 -
               m.m21 * m.m03 * m.m32 -
               m.m31 * m.m02 * m.m23 +
               m.m31 * m.m03 * m.m22;
    inv.m11 =  m.m00 * m.m22 * m.m33 -
               m.m00 * m.m23 * m.m32 -
               m.m20 * m.m02 * m.m33 +
               m.m20 * m.m03 * m.m32 +
               m.m30 * m.m02 * m.m23 -
               m.m30 * m.m03 * m.m22;
    inv.m21 = -m.m00 * m.m21 * m.m33 +
               m.m00 * m.m23 * m.m31 +
               m.m20 * m.m01 * m.m33 -
               m.m20 * m.m03 * m.m31 -
               m.m30 * m.m01 * m.m23 +
               m.m30 * m.m03 * m.m21;
    inv.m31 =  m.m00 * m.m21 * m.m32 -
               m.m00 * m.m22 * m.m31 -
               m.m20 * m.m01 * m.m32 +
               m.m20 * m.m02 * m.m31 +
               m.m30 * m.m01 * m.m22 -
               m.m30 * m.m02 * m.m21;
    inv.m02 =  m.m01 * m.m12 * m.m33 -
               m.m01 * m.m13 * m.m32 -
               m.m11 * m.m02 * m.m33 +
               m.m11 * m.m03 * m.m32 +
               m.m31 * m.m02 * m.m13 -
               m.m31 * m.m03 * m.m12;
    inv.m12 = -m.m00 * m.m12 * m.m33 +
               m.m00 * m.m13 * m.m32 +
               m.m10 * m.m02 * m.m33 -
               m.m10 * m.m03 * m.m32 -
               m.m30 * m.m02 * m.m13 +
               m.m30 * m.m03 * m.m12;
    inv.m22 =  m.m00 * m.m11 * m.m33 -
               m.m00 * m.m13 * m.m31 -
               m.m10 * m.m01 * m.m33 +
               m.m10 * m.m03 * m.m31 +
               m.m30 * m.m01 * m.m13 -
               m.m30 * m.m03 * m.m11;
    inv.m32 = -m.m00 * m.m11 * m.m32 +
               m.m00 * m.m12 * m.m31 +
               m.m10 * m.m01 * m.m32 -
               m.m10 * m.m02 * m.m31 -
               m.m30 * m.m01 * m.m12 +
               m.m30 * m.m02 * m.m11;
    inv.m03 = -m.m01 * m.m12 * m.m23 +
               m.m01 * m.m13 * m.m22 +
               m.m11 * m.m02 * m.m23 -
               m.m11 * m.m03 * m.m22 -
               m.m21 * m.m02 * m.m13 +
               m.m21 * m.m03 * m.m12;
    inv.m13 =  m.m00 * m.m12 * m.m23 -
               m.m00 * m.m13 * m.m22 -
               m.m10 * m.m02 * m.m23 +
               m.m10 * m.m03 * m.m22 +
               m.m20 * m.m02 * m.m13 -
               m.m20 * m.m03 * m.m12;
    inv.m23 = -m.m00 * m.m11 * m.m23 +
               m.m00 * m.m13 * m.m21 +
               m.m10 * m.m01 * m.m23 -
               m.m10 * m.m03 * m.m21 -
               m.m20 * m.m01 * m.m13 +
               m.m20 * m.m03 * m.m11;
    inv.m33 =  m.m00 * m.m11 * m.m22 -
               m.m00 * m.m12 * m.m21 -
               m.m10 * m.m01 * m.m22 +
               m.m10 * m.m02 * m.m21 +
               m.m20 * m.m01 * m.m12 -
               m.m20 * m.m02 * m.m11;

    float det = m.m00 * inv.m00 + m.m01 * inv.m10 + m.m02 * inv.m20 + m.m03 * inv.m30;
    det = 1 / det;

    for (int i = 0; i < 16; ++i) {
        ((float *) &inv)[i] *= det;
    }

    return inv;

    //212 multiplies
    //83 adds
    //1 divide
}

//this imlementation is adapted from GLM
inline Mat4 inverse_transpose(Mat4 m) {
    //sub-factors
    float f00 = m.m22 * m.m33 - m.m32 * m.m23;
    float f01 = m.m21 * m.m33 - m.m31 * m.m23;
    float f02 = m.m21 * m.m32 - m.m31 * m.m22;
    float f03 = m.m20 * m.m33 - m.m30 * m.m23;
    float f04 = m.m20 * m.m32 - m.m30 * m.m22;
    float f05 = m.m20 * m.m31 - m.m30 * m.m21;
    float f06 = m.m12 * m.m33 - m.m32 * m.m13;
    float f07 = m.m11 * m.m33 - m.m31 * m.m13;
    float f08 = m.m11 * m.m32 - m.m31 * m.m12;
    float f09 = m.m10 * m.m33 - m.m30 * m.m13;
    float f10 = m.m10 * m.m32 - m.m30 * m.m12;
    float f11 = m.m11 * m.m33 - m.m31 * m.m13;
    float f12 = m.m10 * m.m31 - m.m30 * m.m11;
    float f13 = m.m12 * m.m23 - m.m22 * m.m13;
    float f14 = m.m11 * m.m23 - m.m21 * m.m13;
    float f15 = m.m11 * m.m22 - m.m21 * m.m12;
    float f16 = m.m10 * m.m23 - m.m20 * m.m13;
    float f17 = m.m10 * m.m22 - m.m20 * m.m12;
    float f18 = m.m10 * m.m21 - m.m20 * m.m11;

    Mat4 inv;
    inv.m00 = + (m.m11 * f00 - m.m12 * f01 + m.m13 * f02);
    inv.m01 = - (m.m10 * f00 - m.m12 * f03 + m.m13 * f04);
    inv.m02 = + (m.m10 * f01 - m.m11 * f03 + m.m13 * f05);
    inv.m03 = - (m.m10 * f02 - m.m11 * f04 + m.m12 * f05);
    inv.m10 = - (m.m01 * f00 - m.m02 * f01 + m.m03 * f02);
    inv.m11 = + (m.m00 * f00 - m.m02 * f03 + m.m03 * f04);
    inv.m12 = - (m.m00 * f01 - m.m01 * f03 + m.m03 * f05);
    inv.m13 = + (m.m00 * f02 - m.m01 * f04 + m.m02 * f05);
    inv.m20 = + (m.m01 * f06 - m.m02 * f07 + m.m03 * f08);
    inv.m21 = - (m.m00 * f06 - m.m02 * f09 + m.m03 * f10);
    inv.m22 = + (m.m00 * f11 - m.m01 * f09 + m.m03 * f12);
    inv.m23 = - (m.m00 * f08 - m.m01 * f10 + m.m02 * f12);
    inv.m30 = - (m.m01 * f13 - m.m02 * f14 + m.m03 * f15);
    inv.m31 = + (m.m00 * f13 - m.m02 * f16 + m.m03 * f17);
    inv.m32 = - (m.m00 * f14 - m.m01 * f16 + m.m03 * f18);
    inv.m33 = + (m.m00 * f15 - m.m01 * f17 + m.m02 * f18);

    float det =
        + m.m00 * inv.m00
        + m.m01 * inv.m01
        + m.m02 * inv.m02
        + m.m03 * inv.m03;

    det = 1 / det;

    for (int i = 0; i < 16; ++i) {
        ((float *) &inv)[i] *= det;
    }

    return inv;

    //104 multiplies
    //61 adds
    //1 divide
}

inline Mat4 ortho(float left, float right, float bottom, float top, float near, float far) {
    Mat4 ret = {};
    ret.m00 =  2 / (right - left);
    ret.m11 =  2 / (top - bottom);
    ret.m22 = -2 / (far - near);
    ret.m33 =  1;
    ret.m30 = - (right + left) / (right - left);
    ret.m31 = - (top + bottom) / (top - bottom);
    ret.m32 = - (far + near) / (far - near);
    return ret;
}

