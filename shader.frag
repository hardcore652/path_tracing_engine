// TODO:
//	  Fix reflection artifacts (probably use another functions):
//    	Fix problems with ellipsoid reflections
//    	Fix torus artifacts
//	  New primitives
//    Low quality version of engine - use reflectiveness parameter without mixing with diffuse reflection
//    Fix errors on nVidia GPUs
//    Optimize the engine

#version 460
#extension GL_EXT_gpu_shader4: enable

uniform vec2 u_resolution;
uniform sampler2D u_tex;
uniform float u_time;
uniform vec2 u_mouse;
uniform vec3 u_cameraPos;
uniform vec2 u_seed;
uniform vec2 u_seed2;
uniform float u_samplePart;
uniform sampler2D u_skybox;
uniform float u_skyBrightness;
uniform float u_cameraAperture;
uniform float u_focalDist;

vec3 sunDir = normalize(vec3(0.18, -0.62, -0.82));
vec3 sunCol = vec3(1.0, 0.98, 0.5);
float planeZ = 0.0;
vec4 planeColor = vec4(0.6, 0.9, 0.6, 0.0001);
vec3 planeNormal = vec3(0.0, 0.0, -1.0);
const int spheresNum = 23;
const int objects2Num = 25;

const float antiAliasing = 0.5;
const bool antiAliasing_enabled = true;
const bool dof_enabled = true;
const bool skybox_enabled = true;
const float MAXDISTANCE = 500.0;
const int MAXREFLECTIONS = 8;
const int samplesCount = 10;
const float FOV = 1.0;
const float sunBrightness = 55.0;

float focalDist = u_focalDist;

mat2 rot(in float a)
{
	float s = sin(a);
	float c = cos(a);
	return mat2(c, -s, s, c);
}

uvec4 R_STATE;
uint TausStep(uint z, int S1, int S2, int S3, uint M) {
	uint b = (((z << S1) ^ z) >> S2);
	return (((z & M) << S3) ^ b);	
}
uint LCGStep(uint z, uint A, uint C) { return (A * z + C); }
vec2 hash22(vec2 p) {
	p += u_seed.x;
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
	p3 += dot(p3, p3.yzx+33.33);
	return fract((p3.xx+p3.yz)*p3.zy);
}
float random() {
	R_STATE.x = TausStep(R_STATE.x, 13, 19, 12, uint(4294967294));
	R_STATE.y = TausStep(R_STATE.y, 2, 25, 4, uint(4294967288));
	R_STATE.z = TausStep(R_STATE.z, 3, 11, 17, uint(4294967280));
	R_STATE.w = LCGStep(R_STATE.w, uint(1664525), uint(1013904223));
	return 2.3283064365387e-10 * float((R_STATE.x ^ R_STATE.y ^ R_STATE.z ^ R_STATE.w));
}
void configureRandom(in vec2 fragPos) {
	vec2 uvRes = hash22(fragPos + 1.0) * u_resolution + u_resolution;
	R_STATE.x = uint(u_seed.x + uvRes.x);
	R_STATE.y = uint(u_seed.y + uvRes.x);
	R_STATE.z = uint(u_seed2.x + uvRes.y);
	R_STATE.w = uint(u_seed2.y + uvRes.y);
}
vec3 randomOnSphere() {
	vec3 rand = vec3(random(), random(), random());
	float theta = rand.x * 2.0 * 3.14159265;
	float v = rand.y;
	float phi = acos(2.0 * v - 1.0);
	float r = pow(rand.z, 1.0 / 3.0);
	float x = r * sin(phi) * cos(theta);
	float y = r * sin(phi) * sin(theta);
	float z = r * cos(phi);
	return vec3(x, y, z);
}

float capIntersect(in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float ra) {
	vec3 ba = pb - pa;
	vec3 oa = ro - pa;
	float baba = dot(ba, ba);
	float bard = dot(ba, rd);
	float baoa = dot(ba, oa);
	float rdoa = dot(rd, oa);
	float oaoa = dot(oa, oa);
	float a = baba - bard * bard;
	float b = baba * rdoa - baoa * bard;
	float c = baba * oaoa - baoa * baoa - ra * ra * baba;
	float h = b * b - a * c;
	if (h >= 0.0) {
		float t = (-b - sqrt(h)) / a;
		float y = baoa + t * bard;
		if (y > 0.0 && y < baba) return t;
		vec3 oc = (y <= 0.0) ? oa:ro - pb;
		b = dot(rd, oc);
		c = dot(oc, oc) - ra * ra;
		h = b * b - c;
		if (h > 0.0) return -b - sqrt(h);
	}
	return -1.0;
}
vec3 capNormal(in vec3 pos, in vec3 a, in vec3 b, in float r) {
	vec3  ba = b - a;
	vec3  pa = pos - a;
	float h  = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
	return (pa - h * ba) / r;
}

vec2 hexPrismInters(in vec3  ro, in vec3  rd, in float ra, in float he, out vec3 n) {
	const float ks3 = 0.866025;
	const vec3 n1 = vec3( 1.0, 0.0, 0.0);
	const vec3 n2 = vec3( 0.5, 0.0, ks3);
	const vec3 n3 = vec3(-0.5, 0.0, ks3);
	const vec3 n4 = vec3( 0.0, 1.0, 0.0);
	vec3 t1 = vec3((vec2(ra, -ra) - dot(ro,n1)) / dot(rd,n1), 1.0);
	vec3 t2 = vec3((vec2(ra, -ra) - dot(ro,n2)) / dot(rd,n2), 1.0);
	vec3 t3 = vec3((vec2(ra, -ra) - dot(ro,n3)) / dot(rd,n3), 1.0);
	vec3 t4 = vec3((vec2(he, -he) - dot(ro,n4)) / dot(rd,n4), 1.0);
	if (t1.y < t1.x) t1 = vec3(t1.yx, -1.0);
	if (t2.y < t2.x) t2 = vec3(t2.yx, -1.0);
	if (t3.y < t3.x) t3 = vec3(t3.yx, -1.0);
	if (t4.y < t4.x) t4 = vec3(t4.yx, -1.0);
	vec4 tN = vec4(t1.x, t1.z * n1);
	if (t2.x > tN.x) tN = vec4(t2.x, t2.z * n2);
	if (t3.x > tN.x) tN = vec4(t3.x, t3.z * n3);
	if (t4.x > tN.x) tN = vec4(t4.x, t4.z * n4);
	float tF = min(min(t1.y, t2.y), min(t3.y, t4.y));
	if (tN.x > tF || tF < 0.0) return vec2(-1.0);
	n = tN.yzw;
	return vec2(tN.x, tF);
}

vec2 eliIntersect(in vec3 ro, in vec3 rd, in vec3 ra) {
    vec3 ocn = ro / ra;
    vec3 rdn = rd / ra;
    float a = dot(rdn, rdn);
    float b = dot(ocn, rdn);
    float c = dot(ocn, ocn);
    float h = b * b - a * (c - 1.0);
    if (h < 0.0) return vec2(-1.0); // no intersection
    h = sqrt(h);
    return vec2(-b - h, -b + h) / a;
}

float torIntersect(in vec3 ro, in vec3 rd, in vec2 tor) {
    float po = 1.0;
    float Ra2 = tor.x * tor.x;
    float ra2 = tor.y * tor.y;
    float m = dot(ro, ro);
    float n = dot(ro, rd);
    float k = (m + Ra2 - ra2) / 2.0;
    float k3 = n;
    float k2 = n * n - Ra2 * dot(rd.xy, rd.xy) + k;
    float k1 = n * k - Ra2 * dot(rd.xy, ro.xy);
    float k0 = k * k - Ra2 * dot(ro.xy, ro.xy);
    if (abs(k3 * (k3 * k3 - k2) + k1) < 0.01) {
        po = -1.0;
        float tmp = k1; k1 = k3; k3 = tmp;
        k0 = 1.0 / k0;
        k1 = k1 * k0;
        k2 = k2 * k0;
        k3 = k3 * k0;
    }
    float c2 = k2 * 2.0 - 3.0 * k3 * k3;
    float c1 = k3 * (k3 * k3 - k2) + k1;
    float c0 = k3 * (k3 * (c2 + 2.0 * k2) - 8.0 * k1) + 4.0 * k0;
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;
    float Q = c2 * c2 + c0;
    float R = c2 * c2 * c2 - 3.0 * c2 * c0 + c1 * c1;
    float h = R * R - Q * Q * Q;
    if (h >= 0.0) {
        h = sqrt(h);
        float v = sign(R + h) * pow(abs(R + h), 1.0 / 3.0); // cube root
        float u = sign(R - h) * pow(abs(R - h), 1.0 / 3.0); // cube root
        vec2 s = vec2((v + u) + 4.0 * c2, (v - u) * sqrt(3.0));
        float y = sqrt(0.5 * (length(s) + s.x));
        float x = 0.5 * s.y / y;
        float r = 2.0 * c1 / (x * x + y * y);
        float t1 =  x - r - k3; t1 = (po < 0.0) ? 2.0 / t1:t1;
        float t2 = -x - r - k3; t2 = (po < 0.0) ? 2.0 / t2:t2;
        float t = 1e20;
        if (t1 > 0.0) t = t1;
        if (t2 > 0.0) t = min(t, t2);
        return t;
    }
    float sQ = sqrt(Q);
    float w = sQ * cos(acos(-R / (sQ * Q)) / 3.0);
    float d2 = -(w + c2); if(d2 < 0.0) return -1.0;
    float d1 = sqrt(d2);
    float h1 = sqrt(w - 2.0 * c2 + c1 / d1);
    float h2 = sqrt(w - 2.0 * c2 - c1 / d1);
    float t1 = -d1 - h1 - k3; t1 = (po < 0.0) ? 2.0 / t1:t1;
    float t2 = -d1 + h1 - k3; t2 = (po < 0.0) ? 2.0 / t2:t2;
    float t3 =  d1 - h2 - k3; t3 = (po < 0.0) ? 2.0 / t3:t3;
    float t4 =  d1 + h2 - k3; t4 = (po < 0.0) ? 2.0 / t4:t4;
    float t = 1e20;
    if (t1 > 0.0) t = t1;
    if (t2 > 0.0) t = min(t, t2);
    if (t3 > 0.0) t = min(t, t3);
    if (t4 > 0.0) t = min(t, t4);
    return t;
}
vec3 torNormal(in vec3 pos, vec2 tor) {
    return normalize(pos * (dot(pos, pos) - tor.y * tor.y - tor.x * tor.x * vec3(1.0, 1.0, -1.0)));
}

float roundedboxIntersect(in vec3 ro, in vec3 rd, in vec3 size, in float rad) {
	vec3 m = 1.0 / rd;
	vec3 n = m * ro;
	vec3 k = abs(m) * (size + rad);
	vec3 t1 = -n - k;
	vec3 t2 = -n + k;
	float tN = max(max(t1.x, t1.y ), t1.z);
	float tF = min(min(t2.x, t2.y ), t2.z);
	if(tN > tF || tF < 0.0) return -1.0;
	float t = tN;
	vec3 pos = ro+t*rd;
	vec3 s = sign(pos);
	ro *= s; rd *= s;
	pos *= s; pos -= size;
	pos = max(pos.xyz, pos.yzx);
	if(min(min(pos.x,pos.y),pos.z) < 0.0) return t;
	vec3 oc = ro - size;
	vec3 dd = rd*rd;
	vec3 oo = oc*oc;
	vec3 od = oc*rd;
	float ra2 = rad*rad;
	t = 1e20;        
	{
		float b = od.x + od.y + od.z;
		float c = oo.x + oo.y + oo.z - ra2;
		float h = b * b - c;
		if (h > 0.0) t = -b - sqrt(h); }
    {
		float a = dd.y + dd.z;
		float b = od.y + od.z;
		float c = oo.y + oo.z - ra2;
		float h = b * b - a * c;
		if (h > 0.0) {
			h = (-b - sqrt(h)) / a;
			if (h > 0.0 && h < t && abs(ro.x + rd.x * h) < size.x) t = h;
		} }
	{
		float a = dd.z + dd.x;
		float b = od.z + od.x;
		float c = oo.z + oo.x - ra2;
		float h = b * b - a * c;
		if (h > 0.0) {
			h = (-b - sqrt(h)) / a;
			if (h > 0.0 && h < t && abs(ro.y + rd.y * h) < size.y) t = h;
		} }
	{
		float a = dd.x + dd.y;
		float b = od.x + od.y;
		float c = oo.x + oo.y - ra2;
		float h = b * b - a * c;
		if (h > 0.0) {
			h = (-b - sqrt(h)) / a;
			if (h > 0.0 && h < t && abs(ro.z + rd.z * h) < size.z) t = h;
		} }
	if (t > 1e19) t = -1.0;  
	return t;
}
vec3 roundedboxNormal(in vec3 pos, in vec3 siz, in float rad) {
	return sign(pos) * normalize(max(abs(pos) - siz, 0.0));
}

vec2 sphIntersect(in vec3 ro, in vec3 rd, in float ra)
{
	float b = dot(ro, rd);
	float c = dot(ro, ro) - ra * ra;
	float h = b * b - c;
	if (h < 0.0) return vec2(-1.0); // no intersection
	h = sqrt(h);
	return vec2(-b - h, -b + h);
}
vec2 boxIntersection(in vec3 ro, in vec3 rd, vec3 boxSize, out vec3 outNormal) 
{
	vec3 m = 1.0 / rd; // can precompute if traversing a set of aligned boxes
	vec3 n = m * ro;   // can precompute if traversing a set of aligned boxes
	vec3 k = abs(m) * boxSize;
	vec3 t1 = -n - k;
	vec3 t2 = -n + k;
	float tN = max(max(t1.x, t1.y), t1.z);
	float tF = min(min(t2.x, t2.y), t2.z);
	if (tN > tF || tF < 0.0) return vec2(-1.0); // no intersection
	outNormal = -sign(rd) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
	return vec2(tN, tF);
}
float plaIntersect(in vec3 ro, in vec3 rd, in vec4 p) { return -(dot(ro, p.xyz) + p.w) / dot(rd, p.xyz); }

vec3 getSky(in vec3 rayDir) 
{
	vec2 fragPos = vec2(atan(rayDir.x, rayDir.y), asin(rayDir.z) * 2.0);
	fragPos /= 3.14159265;
	fragPos *= 0.5;
	fragPos += 0.5;
	vec3 sky = texture(u_skybox, fragPos).rgb;

	vec3 sun = sunCol;
	sun *= max(0.0, pow(dot(rayDir, sunDir), 2048.0));
	sun *= sunBrightness;
	sun += sunCol * max(0.0, pow(dot(rayDir, sunDir), 30.0)) / 240.0;
	vec3 color = sun + sky * 0.01;
	return color / u_skyBrightness;
}

vec2 simpleRayCast(in vec3 rayOrig, in vec3 rayDir, in vec4 spheres[spheresNum][3], in vec3 boxes[objects2Num][5]) {
	vec2 minInters = vec2(MAXDISTANCE);
	vec2 inters;
	for (int i = 0; i < spheres.length(); i++) {
		inters = sphIntersect(rayOrig - spheres[i][0].xyz, rayDir, spheres[i][0].w);
		if (inters.x > 0.0 && inters.x < minInters.x) minInters = inters; }
	for (int i = 0; i < boxes.length(); i++) {
		vec3 ro = rayOrig - boxes[i][0];
		vec3 rd = rayDir;
		if (boxes[i][4].z != 0.0) {
			ro.xy *= rot(boxes[i][4].z);
			rd.xy *= rot(boxes[i][4].z); }
		if (boxes[i][4].x != 0.0) {
			ro.yz *= rot(boxes[i][4].x);
			rd.yz *= rot(boxes[i][4].x); }
		if (boxes[i][4].y != 0.0) {
			ro.xz *= rot(boxes[i][4].y);
			rd.xz *= rot(boxes[i][4].y); }
		vec3 n;
		if (boxes[i][3].z == -1.0) inters = eliIntersect(ro, rd, boxes[i][1]);
		else if (boxes[i][3].z == -2.0) inters = hexPrismInters(ro, rd, boxes[i][1].x, boxes[i][1].z, n);
		else if (boxes[i][3].z == -3.0) inters = vec2(torIntersect(ro, rd, boxes[i][1].xz));
		else if (boxes[i][3].z == -4.0) inters = vec2(capIntersect(rayOrig, rayDir, boxes[i][0], boxes[i][4], boxes[i][1].x));
		else inters = boxIntersection(ro, rd, boxes[i][1], n);
		if (inters.x > 0.0 && inters.x < minInters.x) minInters = inters; }
	inters = vec2(plaIntersect(rayOrig - vec3(0.0, 0.0, planeZ), rayDir, vec4(planeNormal, 1.0)));
	if (inters.x > 0.0 && inters.x < minInters.x) minInters = inters;
	if (minInters == vec2(MAXDISTANCE)) minInters = vec2(20.0);
	return minInters;
}

vec4 rayCast(inout vec3 rayOrig, inout vec3 rayDir, in vec4 spheres[spheresNum][3], in vec3 boxes[objects2Num][5]) {
	vec2 minInters = vec2(MAXDISTANCE);
	vec2 inters;
	vec3 norm;
	vec3 startPos;
	vec4 color;
	vec4 additMat = vec4(0.0);
	
	for (int i = 0; i < spheres.length(); i++) {
		startPos = rayOrig - spheres[i][0].xyz;
		inters = sphIntersect(startPos, rayDir, spheres[i][0].w);
		if (inters.x > 0.0 && inters.x < minInters.x) {
			minInters = inters;
			norm = normalize(startPos + rayDir * inters.x);
			color = spheres[i][1];
			additMat = spheres[i][2];
		}
	}
	for (int i = 0; i < boxes.length(); i++) {
		vec3 tempNorm;
		if (boxes[i][3].z != -4.0) {
			vec3 ro = rayOrig - boxes[i][0];
			vec3 rd = rayDir;
			if (boxes[i][4].z != 0.0) {
				ro.xy *= rot(boxes[i][4].z);
				rd.xy *= rot(boxes[i][4].z); }
			if (boxes[i][4].x != 0.0) {
				ro.yz *= rot(boxes[i][4].x);
				rd.yz *= rot(boxes[i][4].x); }
			if (boxes[i][4].y != 0.0) {
				ro.xz *= rot(boxes[i][4].y);
				rd.xz *= rot(boxes[i][4].y); }

			if (boxes[i][3].z == 0.0) inters = boxIntersection(ro, rd, boxes[i][1], tempNorm);
			else if (boxes[i][3].z == -1.0) {
				inters = eliIntersect(ro, rd, boxes[i][1]);
				tempNorm = normalize((ro + rd * inters.x) / boxes[i][1]);
			}
			else if (boxes[i][3].z == -2.0) inters = hexPrismInters(ro, rd, boxes[i][1].x, boxes[i][1].z, tempNorm);
			else if (boxes[i][3].z == -3.0) {
				inters = vec2(torIntersect(ro, rd, boxes[i][1].xz));
				tempNorm = torNormal(ro + rd * inters.x, boxes[i][1].xz);
			}
			else if (boxes[i][3].z > 0.0) {
				inters = vec2(roundedboxIntersect(ro, rd, boxes[i][1]-boxes[i][3].z/2.0, boxes[i][3].z));
				tempNorm = roundedboxNormal(ro + rd * inters.x, boxes[i][1]-boxes[i][3].z/2.0, boxes[i][3].z);
			}
			if (boxes[i][4].y != 0.0) tempNorm.xz *= rot(-boxes[i][4].y);
			if (boxes[i][4].x != 0.0) tempNorm.yz *= rot(-boxes[i][4].x);
			if (boxes[i][4].z != 0.0) tempNorm.xy *= rot(-boxes[i][4].z);
		}
		else {
			inters = vec2(capIntersect(rayOrig, rayDir, boxes[i][0], boxes[i][4], boxes[i][1].x));
			tempNorm = capNormal(rayOrig + rayDir * inters.x, boxes[i][0], boxes[i][4], boxes[i][1].x);
		}
		
		if (inters.x > 0.0 && inters.x < minInters.x) {
			minInters = inters;
			norm = tempNorm;
			color = vec4(boxes[i][2], boxes[i][3].x);
			additMat = vec4(0.0);
			additMat.y = boxes[i][3].y;
		}
	}	
	startPos = rayOrig - vec3(0.0, 0.0, planeZ);
	inters = vec2(plaIntersect(startPos, rayDir, vec4(planeNormal, 1.0)));
	if (inters.x > 0.0 && inters.x < minInters.x) {
		additMat = vec4(0.0);
		minInters = inters;
		norm = planeNormal;
		color = planeColor;
	}
	if (minInters.x == MAXDISTANCE) {
		if (skybox_enabled) return vec4(getSky(rayDir), -1.0);
		else return vec4(vec3(0.003, 0.006, 0.008) / u_skyBrightness, -1.0); }
	if (color.a == -1.0) return color;
	if (color.a <= -2.0) {
		if (random() - 3.0 - color.a > 0.0) color.a = 1.0;
		else color.a = 0.0;
	}
	vec3 reflected = reflect(rayDir, norm);

	if (additMat.x > 0.0) {
		float fresnel = 1.0 - abs(dot(-rayDir, norm));
		if(random() - 0.1 - color.a < fresnel * fresnel) {
			rayDir = reflected;
			return color;
		}
		rayOrig += rayDir * (minInters.y + 0.001);
		rayDir = refract(rayDir, norm, additMat.x);
		return color;
	}
	if (additMat.y > 0.0) {
		if (random() < additMat.y) {
			rayOrig += rayDir * (minInters.y + 0.001);
			return color;
		}
	}

	vec3 intersPos = rayOrig + rayDir * inters.x;
	vec3 rand = randomOnSphere();
	vec3 diffuse = normalize(rand * dot(rand, norm));
	rayOrig += rayDir * (minInters.x - 0.01);
	rayDir = mix(diffuse, reflected, color.a);
	return color;
}

vec3 traceRay(in vec3 rayOrig, in vec3 rayDir, in vec4 spheres[spheresNum][3], in vec3 objects2[objects2Num][5]) {
	vec3 color = vec3(1.0);
	for (int i = 0; i < MAXREFLECTIONS; i++)
	{
		vec4 reflectionCol = rayCast(rayOrig, rayDir, spheres, objects2);
		color *= reflectionCol.rgb;
		if(reflectionCol.a == -1.0) return color;
	}
	return vec3(0.0);
	

}

vec2 randomPointOnCircle(float R) {
	float r = R * sqrt(random());
	float theta = random() * 2 * 3.14159265;
	return vec2(r * cos(theta), r * sin(theta));
}

vec2 anti_aliasing(in vec2 uv) {
	vec2 jitter = vec2(random(), random());
	jitter = jitter * 2.0 - 1.0;
	return uv + jitter / u_resolution * antiAliasing;
}

vec3 depth_of_field(in vec2 fragPos, inout vec3 rayDir, in vec3 rayOrig, in float focalDist) {
	vec3 point = vec3(0.0, randomPointOnCircle(u_cameraAperture * focalDist / 10.0));
	rayDir = normalize(vec3(focalDist * FOV, fragPos * focalDist - point.yz));
	point.zx *= rot(-u_mouse.y);
	point.xy *= rot(u_mouse.x);
	rayOrig += point;
	return rayOrig;
}

vec3 colCorr(in vec3 color) {
	float white = 20.0;
	color *= white * 16.0;
	vec3 sampleCol = texture(u_tex, gl_FragCoord.xy / u_resolution).rgb;
	return mix(sampleCol, (color * (1.0 + color / white / white)) / (1.0 + color) * 1.2, u_samplePart);
}

vec4 spheres[spheresNum][3];
vec3 objects2[objects2Num][5];

void main() {

	// Spheres
	// [0] - pos and radius
	// [1] - color and roughness
	// [2] - other materials: x - refractiveness, y - alpha, z - nothing, w - nothing

	spheres[0][0] = vec4(0.0, 0.0, 0.0, 1.0);
	spheres[0][1] = vec4(1.0, 0.8, 0.1, 1.0);

	spheres[1][0] = vec4(4.0, -5.0, -0.51, 1.5);
	spheres[1][1] = vec4(1.0, 1.0, 1.0, 0.0);
	spheres[1][2] = vec4(0.5, 0.0, 0.0, 0.0);

	spheres[2][0] = vec4(-8.0, 5.0, -1.0, 2.0);
	spheres[2][1] = vec4(0.8, 0.8, 1.0, 0.0);

	spheres[3][0] = vec4(-4.0, -6.0, 0.0, 1.0);
	spheres[3][1] = vec4(0.8, 0.8, 1.0, 0.55);

	spheres[4][0] = vec4(-5.0, 3.0, -3.3, 0.7);
	spheres[4][1] = vec4(2.0, 0.004, 0.05, -1.0);

	spheres[5][0] = vec4(-5.0, -2.0, 0.1, 0.9);
	spheres[5][1] = vec4(0.3, 0.3, 0.3, -1.0);

	for (int i = 0; i < 5; i++) {
		spheres[i+6][0] = vec4(10.0 + i * 1.7, 10.0, 0.2, 0.8);
		spheres[i+6][1] = vec4(0.6, 0.6, 0.6, i / 4.0);

	}
	for (int i = 0; i < 5; i++) {
		spheres[i+11][0] = vec4(10.0 + i * 1.7, 12.0, 0.2, 0.8);
		spheres[i+11][1] = vec4(0.6, 0.6, 0.6, -2.0 - i / 4.0);

	}

	spheres[16][0] = vec4(8.0, -3.0, -3.61, 0.4);
	spheres[16][1] = vec4(0.003, 0.009, 0.5, -1.0);
	spheres[16][1].b *= 50.0;

	spheres[17][0] = vec4(-7.0, -6.0, 0.2, 0.8);
	spheres[17][1] = vec4(0.6, 0.6, 1.0, 0.2);
	spheres[17][2] = vec4(0.5, 0.0, 0.0, 0.0);

	spheres[18][0] = vec4(-34.0, -36.0, -3.0, 2.0);
	spheres[18][1] = vec4(1.0, 1.0, 1.0, 1.0);

	spheres[19][0] = vec4(-26.5, -36.7, -4.275, 0.6);
	spheres[19][1] = vec4(0.4, 1.0, 0.45, -2.5);

	spheres[20][0] = vec4(-27.3, -31.0, -2.0, 1.0);
	spheres[20][1] = vec4(1.0, 1.0, 1.0, 0.0);
	spheres[20][2] = vec4(0.5, 0.0, 0.0, 0.0);

	spheres[21][0] = vec4(-30.0, -29.0, -1.8, 0.8);
	spheres[21][1] = vec4(1.0, 0.8, 0.85, 0.01);

	spheres[22][0] = vec4(-34.5, -30.0, -2.1, 1.1);
	spheres[22][1] = vec4(1.0, 1.0, 1.0, 0.5);
	spheres[22][2] = vec4(0.0, 0.2, 0.0, 0.0);



	// Other objects (default - box)
	// [0] - pos   (for capsule: first point)
	// [1] - size, (for hex prism X = radius, Y = nothing, Z = height,
	//              for torus X = size, Z = width,
	//				for capsule X = radius)
	// [2] - color
	// [3] - materials, x - roughness, y - alpha, z - roundness,
	//       if z = -1.0: ellipsoid,
	//       if z = -2.0: hexagonal prism,
	//       if z = -3.0: torus,
	//		 if z = -4.0: capsule
	// [4] - rotation (for capsule: second point)


	objects2[0][0] = vec3(0.0, 3.0, 0.5001);
	objects2[0][1] = vec3(1.0, 1.3, 1.5);
	objects2[0][2] = vec3(0.1, 0.1, 1.0);
	objects2[0][3] = vec3(0.3, 0.0, 0.0);
	objects2[0][4] = vec3(0.0, 0.0, 0.1);

	objects2[1][0] = vec3(8.0, -3.0, -0.01);
	objects2[1][1] = vec3(1.0, 1.0, 1.0);
	objects2[1][2] = vec3(0.6, 0.6, 0.6);
	objects2[1][3] = vec3(0.7, 0.2, 0.0);
	objects2[1][4] = vec3(0.0, 0.0, 0.2);

	objects2[3][0] = vec3(8.0, -3.0, -1.71);
	objects2[3][1] = vec3(0.7, 0.7, 0.7);
	objects2[3][2] = vec3(0.6, 0.6, 0.6);
	objects2[3][3] = vec3(0.7, 0.2, 0.0);
	objects2[3][4] = vec3(0.0, 0.0, 0.5);

	objects2[4][0] = vec3(8.0, -3.0, -2.81);
	objects2[4][1] = vec3(0.4, 0.4, 0.4);
	objects2[4][2] = vec3(0.6, 0.6, 0.6);
	objects2[4][3] = vec3(0.7, 0.2, 0.0);
	objects2[4][4] = vec3(0.0, 0.0, 0.7);

	objects2[2][0] = vec3(0.0, 3.0, -1.3);
	objects2[2][1] = vec3(0.4, 0.4, 0.1);
	objects2[2][2] = vec3(0.008, 0.3, 0.008);
	objects2[2][3] = vec3(-1.0, 0.0, 0.0);

	objects2[5][0] = vec3(8.0, 17.0, -6.00);	
	objects2[5][1] = vec3(2.0, 2.0, 0.15);
	objects2[5][2] = vec3(0.006, 3.0, 0.006);
	objects2[5][3] = vec3(-1.0, 0.0, 0.0);
	objects2[5][4] = vec3(0.0, 0.6, 0.9);
													
	objects2[6][0] = vec3(21.0, 6.0, -6.00);			// puple light
	objects2[6][1] = vec3(2.0, 2.0, 0.15);
	objects2[6][2] = vec3(3.0, 0.006, 0.1);
	objects2[6][3] = vec3(-1.0, 0.0, 0.0);
	objects2[6][4] = vec3(0.0, -0.6, 0.9);

	objects2[7][0] = vec3(30.0, 30.0, -1.8);			// green light
	objects2[7][1] = vec3(2.0, 2.0, 2.0);
	objects2[7][2] = vec3(1.0, 1.0, 1.0);
	objects2[7][3] = vec3(0.005, 0.0, 0.0);
	objects2[7][4] = vec3(0.7853981625, 0.0, 2.0);

	objects2[8][0] = vec3(-38.0, -30.0, -5.0);
	objects2[8][1] = vec3(1.0, 1.0, 0.0001);
	objects2[8][2] = vec3(8.0, 8.0, 8.0);
	objects2[8][3] = vec3(-1.0, 0.0, 0.0);
	objects2[8][4] = vec3(0.0, 3.14159265 / 2.0, 0.0);

														// cornell box
	objects2[9][0] = vec3(-30.0, -30.0, 0.0);				// floor	
	objects2[9][1] = vec3(8.0, 8.0, 1.0);
	objects2[9][2] = vec3(1.0, 1.0, 1.0);
	objects2[9][3] = vec3(0.0000001, 0.0, 0.0);

	objects2[10][0] = vec3(-38.0, -30.0, -1.0);				// walls
	objects2[10][1] = vec3(0.0, 8.0, 8.0);
	objects2[10][2] = vec3(1.0, 0.87, 0.45);
	objects2[10][3] = vec3(0.0000001, 0.0, 0.0);

	objects2[11][0] = vec3(-22.0, -30.0, -1.0);
	objects2[11][1] = vec3(0.0, 8.0, 8.0);
	objects2[11][2] = vec3(0.4, 0.6, 1.0);
	objects2[11][3] = vec3(0.0000001, 0.0, 0.0);

	objects2[12][0] = vec3(-30.0, -38.0, -1.0);
	objects2[12][1] = vec3(8.0, 0.0, 8.0);
	objects2[12][2] = vec3(1.0, 1.0, 1.0);
	objects2[12][3] = vec3(0.0000001, 0.0, 0.0);

	objects2[13][0] = vec3(-30.0, -22.0, -1.0);
	objects2[13][1] = vec3(8.0, 0.0, 8.0);
	objects2[13][2] = vec3(1.0, 1.0, 1.0);
	objects2[13][3] = vec3(0.0000001, 0.0, 0.0);			// end (walls)
														
	objects2[14][0] = vec3(-30.0, -30.0, -9.0);				// top
	objects2[14][1] = vec3(8.0, 8.0, 0.0);
	objects2[14][2] = vec3(1.0, 1.0, 1.0);
	objects2[14][3] = vec3(0.0000001, 0.0, 0.0);

	objects2[15][0] = vec3(-26.5, -36.7, -1.9);				// cube in box
	objects2[15][1] = vec3(1.3, 1.3, 1.75);
	objects2[15][2] = vec3(1.0, 1.0, 1.0);
	objects2[15][3] = vec3(0.0000001, 0.0, 0.0);
														// end (cornell box)

	objects2[16][0] = vec3(-30.0, 30.0, -1.8);			// rounded boxes
	objects2[16][1] = vec3(2.0, 2.0, 2.0);
	objects2[16][2] = vec3(1.0, 1.0, 1.0);
	objects2[16][3] = vec3(0.005, 0.0, 0.5);
	objects2[16][4] = vec3(0.7853981625, 0.0, 2.0);

	objects2[17][0] = vec3(-32.7, 24.0, -0.8);
	objects2[17][1] = vec3(1.0, 1.0, 1.0);
	objects2[17][2] = vec3(1.0, 1.0, 1.0);
	objects2[17][3] = vec3(0.8, 0.3, 0.5);
	objects2[17][4] = vec3(0.7853981625, 0.0, 2.0);		// end (rounded boxes)

	objects2[18][0] = vec3(-34.0, 18.0, 0.0);
	objects2[18][1] = vec3(1.0, 1.0, 1.0);
	objects2[18][2] = vec3(3.0, 3.0, 3.0);
	objects2[18][3] = vec3(-1.0, 0.0, 0.0);

	objects2[19][0] = vec3(-34.0, 0.0, -2.0);			// test ellipsoid
	objects2[19][1] = vec3(2.0, 3.5, 1.7);
	objects2[19][2] = vec3(0.5, 0.5, 0.5);
	objects2[19][3] = vec3(1.0, 0.0, -1.0);
	objects2[19][4] = vec3(0.7853981625, 0.0, 2.0);

	objects2[20][0] = vec3(-34.0, 10.0, -2.0);			// test ellipsoid
	objects2[20][1] = vec3(2.0, 3.0, 2.0);
	objects2[20][2] = vec3(0.5, 0.5, 0.5);
	objects2[20][3] = vec3(0.8, 0.0, -1.0);

	objects2[21][0] = vec3(-34.0, -10.0, -19.0);		// test hex prism
	objects2[21][1] = vec3(2.5, 0.0, 20.0);
	objects2[21][2] = vec3(1.0, 1.0, 1.0);
	objects2[21][3] = vec3(0.001, 0.0, -2.0);
	objects2[21][4] = vec3(3.14159265 / 2.0, 0.0, 0.0);

	objects2[22][0] = vec3(-20.0, 35.0, -1.3);			// test torus
	objects2[22][1] = vec3(2.0, 0.0, 1.0);
	objects2[22][2] = vec3(0.94, 0.9, 0.5);
	objects2[22][3] = vec3(1.0, 0.0, -3.0);
	objects2[22][4] = vec3(1.0, 0.0, 0.0);

	objects2[23][0] = vec3(-12.0, 35.0, -1.3);			// test torus
	objects2[23][1] = vec3(2.0, 0.0, 1.0);
	objects2[23][2] = vec3(0.94, 0.9, 0.5);
	objects2[23][3] = vec3(0.001, 0.0, -3.0);
	objects2[23][4] = vec3(1.0, 0.0, 0.0);

	objects2[24][0] = vec3(2.0, 35.0, -1.3);			// test capsule
	objects2[24][1] = vec3(1.0, 0.0, 0.0);
	objects2[24][2] = vec3(1.0, 0.2, 0.2);
	objects2[24][3] = vec3(-2.5, 0.0, -4.0);
	objects2[24][4] = vec3(-1.0, 35.0, -1.3);


	
	vec2 fragPos = (gl_FragCoord.xy - u_resolution / 2.0) / u_resolution.y; // divide by height to fix diff between screen width and height
	configureRandom(fragPos);
	if (antiAliasing_enabled) fragPos = anti_aliasing(fragPos);

	vec3 rayOrig = u_cameraPos;
	vec3 rayDir;
	if (dof_enabled) {
		if (focalDist == -1.0) {
			vec3 cameraDir = normalize(vec3(1.0, 0.0, 0.0));
			cameraDir.zx *= rot(-u_mouse.y);
			cameraDir.xy *= rot(u_mouse.x);
			focalDist = simpleRayCast(rayOrig, cameraDir, spheres, objects2).x;
		}
		rayOrig = depth_of_field(fragPos, rayDir, rayOrig, focalDist); }
	else rayDir = normalize(vec3(FOV, fragPos));
	rayDir.zx *= rot(-u_mouse.y);
	rayDir.xy *= rot(u_mouse.x);

	vec3 color = vec3(0.0);
	for (int i = 0; i < samplesCount; i++) color += traceRay(rayOrig, rayDir, spheres, objects2);
	color /= samplesCount; 

	gl_FragColor = vec4(colCorr(color), 1.0);
}