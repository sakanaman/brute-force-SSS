#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
struct Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    const double& operator[](const int i) const {return (&x)[i];}
    double& operator[](const int i){return (&x)[i];}
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    double length() const {return sqrt(this->dot(*this));}
};
struct Ray
{
    Vec o, d;
    bool in_medium = false; //始点oが媒質中か、そうでないか。
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    VOL
}; // material types, used in radiance()
struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    {                     // returns distance, 0 if nohit
        Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-8, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
Sphere spheres[] = {
    //Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   //Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), //Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         //Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               //Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         //Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), //Top
    // Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, DIFF),        //Mirr
    Sphere(16.5, Vec(50, 42, 78), Vec(), Vec(1, 1, 1) * .999, VOL),        //Glas
    Sphere(600, Vec(50,681.6-.27,81.6), Vec(12, 12, 12), Vec(), DIFF)     //Lite
};
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

/// utility

Vec Exp(const Vec& v)
{
    return Vec(exp(v.x), exp(v.y), exp(v.z));
}

int Min(int a, int b)
{
    return a < b ? a : b;
}

Vec _calc_normal(const Vec& position, const Vec& centor)
{
    return (position - centor).norm();
}

bool _is_outside(const Vec& wo, const Vec& n)
{
    return wo.dot(n) > 0;//内積がプラスなら外！
}

//正規直交基底
void branchlessONB(const Vec &n, Vec &b1, Vec &b2)
{
    double sign = copysignf(1.0f, n.z);
    const double a = -1.0f / (sign + n.z);
    const double b = n.x * n.y * a;
    b1 = Vec(1.0f + sign * n.x * n.x * a, sign * b, -sign * n.x);
    b2 = Vec(b, sign + n.y * n.y * a, -n.y);
}

Vec local2world(const Vec& local, const Vec& x, const Vec& y, const Vec& z)
{
    return x * local[0] + y * local[1] + z * local[2];
}

// implementation

const Vec sigma_s(0.16, 4.0, 0.24);
const Vec sigma_a(0.008, 0.008, 0.008);
const Vec sigma_t = sigma_s + sigma_a;
double g = 0.1;


Vec Tr(double distance)
{
    Vec inside = sigma_t * distance * (-1); // exp(-sigma_t * dist)
    return Exp(inside);
}

double henyey_greenstein_phase(double theta, double g)
{
    double denomi = 1 + g*g - 2*g*cos(theta);
    return 1/(4 * M_PI) * (1 - g*g) / pow(denomi, 3.0f/2.0f);
}

Vec henyey_greenstein_sample(double g, double u, double v)
{
    double s = 2*u - 1;
    double T = (1 - g*g) / (1 + g * s);
    double cosTheta = 1.0/(2.0 * g) * (1 + g*g - pow(T, 2.0f));
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    double phi = 2 * M_PI * v;

    return {sinTheta * cos(phi),
            sinTheta * sin(phi),
            cosTheta};
}



void Vol_Sample(const Ray& ray, unsigned short *Xi, const double tMax,
               Ray& scatter_Ray, bool& is_scatter, Vec* pdf_throughput, Vec* f_throughput, const int now_channnel)
{
    int channel = now_channnel;

    // distance sampling
    double dist = -log(1 - erand48(Xi))/sigma_t[channel];

    bool sampleMedium = (dist < tMax - 1e-6);

    if(sampleMedium)
    {
        Vec scattering_pos = ray.o + ray.d * dist;
    
        // direction sampling
        Vec loc_axis_x, loc_axis_y;
        branchlessONB(ray.d, loc_axis_x, loc_axis_y);
        Vec local_henyey = henyey_greenstein_sample(g, erand48(Xi), erand48(Xi));
        Vec scattering_dir = local2world(local_henyey, loc_axis_x, loc_axis_y, ray.d);

        //renew scatter_Ray
        scatter_Ray = Ray(scattering_pos, scattering_dir);
        scatter_Ray.in_medium = true;
        is_scatter = true;

        Vec transmit = Tr(dist);
        Vec pdf = sigma_t.mult(transmit);
        Vec f_sss = sigma_s.mult(transmit);

        *pdf_throughput = pdf_throughput->mult(pdf);
        *f_throughput = f_throughput->mult(f_sss);
    }
    else
    {
        Vec transmit = Tr(tMax);
        Vec pdf = transmit;
        Vec f_sss = transmit;
        is_scatter = false;

        *pdf_throughput = pdf_throughput->mult(pdf);
        *f_throughput = f_throughput->mult(f_sss);
    }
}



Vec radience_vol(const Ray &r, int depth, unsigned short *Xi)
{
    Vec pdf_throughput(1.0, 1.0, 1.0);
    Vec f_throughput(1.0, 1.0, 1.0);
    Vec result(0, 0, 0);
    Ray trace_ray = r;
    trace_ray.in_medium = false; //note:最初は媒質外
    double Pr = 1.0;

    int now_channel = Min((int)(erand48(Xi)*3), 2);

    for(int bounds = 0;;bounds++)
    {
        //始点が媒質中
        if(trace_ray.in_medium)
        {
            double t;
            int id = 0;
            if(!intersect(trace_ray, t, id))
            {
                printf("[Caution]: Unhit inside volume");
                break;
            }

            if(spheres[id].refl != VOL)
            {
                printf("%f\n", (trace_ray.o - spheres[6].p).length());
                break;
            }     

            Ray scatterRay(trace_ray.o, trace_ray.d);
            double tMax = t;
            bool is_scatter = false;
            Vol_Sample(trace_ray, Xi, t, scatterRay, is_scatter, 
                        &pdf_throughput, &f_throughput, now_channel);
            
            //散乱する
            if(is_scatter)
            {
                trace_ray = scatterRay;
                trace_ray.in_medium = true;

                continue;
            }

            //散乱しない
            trace_ray = Ray(trace_ray.o + trace_ray.d * t, trace_ray.d);
            trace_ray.in_medium = false;

            continue;
        }
        else //始点が媒質外(要はサーフェスからのレイ)
        {
            double t;
            int id = 0;
            if(!intersect(trace_ray, t, id)) //何も当たんなかった
            {
                break;
            }
            const Sphere &obj = spheres[id];
            // x:衝突点、n:法線、nl:こっち向いてる法線, f:物体の色
            Vec x = trace_ray.o + trace_ray.d * t, n = (x - obj.p).norm(), nl = n.dot(trace_ray.d) < 0 ? n : n * -1, f = obj.c;

            if(obj.e.x > 0)
            {
                double pdf = (pdf_throughput.x + pdf_throughput.y + pdf_throughput.z) * (1.0/3.0);
                if(std::isinf(pdf))
                {
                    break;
                }
                if(std::isinf(f_throughput.x) | std::isinf(f_throughput.x) | std::isinf(f_throughput.x))
                {
                    break;
                }
                result = f_throughput.mult(obj.e) * (1.0/pdf);
                break;
            }

            if(bounds > 5 && erand48(Xi) >= Pr) break;
            pdf_throughput = pdf_throughput * Pr;
            Pr *= 0.98;

            if(obj.refl == DIFF) //拡散反射
            {
                double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
                Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
                Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

                double costheta = abs(d.dot(w));

                trace_ray = Ray(x, d);
                trace_ray.in_medium = false; //もちろん媒質の外
                f_throughput = f_throughput.mult(f);

                continue; //次のターンへ
            }
            else if(obj.refl == VOL) //媒質に当たった
            {
                //もういっちょレイ飛ばす
                Ray involRay(x, trace_ray.d);
                double invol_t;
                int invol_id;
                if(!intersect(involRay, invol_t, invol_id))
                {
                    printf("[Caution]:  [Caution]: Unhit inside volume"); //当たんなかったらバグ
                    break;
                }
                if(spheres[invol_id].refl != VOL)
                {
                    printf("%f\n", (trace_ray.o - spheres[6].p).length());
                    break;
                }     


                double tMax = invol_t;
                Ray scatterRay(x, r.d);
                bool is_scatter = false;
                Vol_Sample(involRay, Xi, tMax, scatterRay, is_scatter, &pdf_throughput, &f_throughput, now_channel);

                if(is_scatter)//散乱する
                {
                    trace_ray = scatterRay;
                    trace_ray.in_medium = true; //もちろん媒質内部

                    continue; //次のターン
                }

                //散乱しない場合
                trace_ray = Ray(involRay.o + involRay.d*invol_t, involRay.d);
                trace_ray.in_medium = false;//次の始点は媒質表面→媒質外とする

                continue; //次のターン
            }
        }
    }
    return result;
}


int main(int argc, char *argv[])
{
    int w = 1024, h = 768, samps = argc == 2 ? atoi(argv[1]) / 4 : 1; // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());        // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
#pragma omp parallel for schedule(dynamic, 1) private(r) // OpenMP
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < w; x++) // Loop cols
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++)       // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radience_vol(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
