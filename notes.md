https://stackoverflow.com/questions/9489736/catmull-rom-curve-with-no-cusps-and-no-self-intersections/23980479#23980479
https://ideone.com/NoEbVM


#include <iostream>
    #include <cmath>
     
    using namespace std;
     
    struct CubicPoly
    {
    	float c0, c1, c2, c3;
     
    	float eval(float t)
    	{
    		float t2 = t*t;
    		float t3 = t2 * t;
    		return c0 + c1*t + c2*t2 + c3*t3;
    	}
    };
     
    /*
     * Compute coefficients for a cubic polynomial
     *   p(s) = c0 + c1*s + c2*s^2 + c3*s^3
     * such that
     *   p(0) = x0, p(1) = x1
     *  and
     *   p'(0) = t0, p'(1) = t1.
     */
    void InitCubicPoly(float x0, float x1, float t0, float t1, CubicPoly &p)
    {
        p.c0 = x0;
        p.c1 = t0;
        p.c2 = -3*x0 + 3*x1 - 2*t0 - t1;
        p.c3 = 2*x0 - 2*x1 + t0 + t1;
    }
     
    // standard Catmull-Rom spline: interpolate between x1 and x2 with previous/following points x1/x4
    // (we don't need this here, but it's for illustration)
    void InitCatmullRom(float x0, float x1, float x2, float x3, CubicPoly &p)
    {
    	// Catmull-Rom with tension 0.5
        InitCubicPoly(x1, x2, 0.5f*(x2-x0), 0.5f*(x3-x1), p);
    }
     
    // compute coefficients for a nonuniform Catmull-Rom spline
    void InitNonuniformCatmullRom(float x0, float x1, float x2, float x3, float dt0, float dt1, float dt2, CubicPoly &p)
    {
        // compute tangents when parameterized in [t1,t2]
        float t1 = (x1 - x0) / dt0 - (x2 - x0) / (dt0 + dt1) + (x2 - x1) / dt1;
        float t2 = (x2 - x1) / dt1 - (x3 - x1) / (dt1 + dt2) + (x3 - x2) / dt2;
     
        // rescale tangents for parametrization in [0,1]
        t1 *= dt1;
        t2 *= dt1;
     
        InitCubicPoly(x1, x2, t1, t2, p);
    }
     
    struct Vec2D
    {
    	Vec2D(float _x, float _y) : x(_x), y(_y) {}
    	float x, y;
    };
     
    float VecDistSquared(const Vec2D& p, const Vec2D& q)
    {
    	float dx = q.x - p.x;
    	float dy = q.y - p.y;
    	return dx*dx + dy*dy;
    }
     
    void InitCentripetalCR(const Vec2D& p0, const Vec2D& p1, const Vec2D& p2, const Vec2D& p3,
    	CubicPoly &px, CubicPoly &py)
    {
        float dt0 = powf(VecDistSquared(p0, p1), 0.25f);
        float dt1 = powf(VecDistSquared(p1, p2), 0.25f);
        float dt2 = powf(VecDistSquared(p2, p3), 0.25f);
     
    	// safety check for repeated points
        if (dt1 < 1e-4f)    dt1 = 1.0f;
        if (dt0 < 1e-4f)    dt0 = dt1;
        if (dt2 < 1e-4f)    dt2 = dt1;
     
    	InitNonuniformCatmullRom(p0.x, p1.x, p2.x, p3.x, dt0, dt1, dt2, px);
    	InitNonuniformCatmullRom(p0.y, p1.y, p2.y, p3.y, dt0, dt1, dt2, py);
    }
     
     
    int main()
    {
    	Vec2D p0(0,0), p1(1,1), p2(1.1,1), p3(2,0);
    	CubicPoly px, py;
    	InitCentripetalCR(p0, p1, p2, p3, px, py);
    	for (int i = 0; i <= 10; ++i)
    		cout << px.eval(0.1f*i) << " " << py.eval(0.1f*i) << endl;
    }



1 1
1.01254 1.00505
1.02346 1.00902
1.03316 1.01189
1.04203 1.01365
1.05046 1.01428
1.05886 1.01377
1.06762 1.0121
1.07713 1.00926
1.08779 1.00523
1.1 1








julia> cmrpoly(u) = [1 u u^2 u^3]
cmrpoly (generic function with 2 methods)

julia> cmrpoints(p1,p2,p3,p4) = Matrix(transpose([p1 p2 p3 p4]))
cmrpoints (generic function with 1 method)

julia> cmrmatrix(τ) = Matrix(transpose(reshape([0 1 0 0 -τ 0 τ 0 2τ τ-3 3-2τ -τ -τ 2-τ τ-2 τ], (4,4))))
cmrmatrix (generic function with 1 method)

julia> cmrpoly(1//7) * cmrmatrix(1//2) * cmrpoints([0,1],[1,5],[2,4],[3,2])
1×2 Array{Rational{Int64},2}:
 8//7  1759//343

