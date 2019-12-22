#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Eigen>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>

using namespace Eigen;
using namespace std;

class Object;

class Camera{
	public:
		map<string, double> bounds;
		Vector3d eye;
		Vector3d look;
		Vector3d up;
		map<string, double> nearFar;
		map<string, double> res;
		Vector3d w;
		Vector3d u;
		Vector3d v;
	Camera()=default;
	void makeCoords(){
		w=eye-look;
		w.normalize();
		u = up.cross(w);
		u.normalize();
		v = w.cross(u);
	}
	
};

class Ray{
	public:
		Vector3d L= Vector3d(-1000,-1000,-1000);
		Vector3d D= Vector3d(-1000,-1000,-1000);
	Ray()=default;
	Ray(Vector3d Lv, Vector3d Dv){
		L = Lv;
		D = Dv;
		D.normalize();
	}
};

class mtl{
	public:
		string name;
		double Ns;
		Vector3d Ks;
		Vector3d Ka;
		Vector3d Kd;
		Vector3d Ke;
		Vector3d Kr;
		double ni;
		double di;
		double illum;
		Vector3d Tr;
	mtl()=default;
	mtl(double Nss, Vector3d Kss, Vector3d Kas, Vector3d Kds, Vector3d Krs,double nis){
		Ns = Nss;
		Ks = Kss;
		Ka = Kas;
		Kd = Kds;
		ni = nis;
		Kr = Krs;
		Vector3d notOpaque(1,1,1); //make sure this is correct later
		if(ni==0){
			Tr = notOpaque;
		}
		else{
			Tr = notOpaque - Kr; 
		}
	}
};

struct HitCatalog {
	double t = numeric_limits<double>::max();
	Vector3d N;
	Vector3d point;
	bool hit = false;
	Object* best;
	mtl mat;
};

class Object{
	public:
		Object()=default;
		virtual bool rayTest(Ray pixel, HitCatalog &h1)=0;
		virtual Vector3d createN(HitCatalog h1)=0;
		virtual Ray refract_exit(Vector3d D, vector<Object*> arrayOfObjects, HitCatalog h1)=0;
};
	
void ray_find(Ray ray, vector<Object*> arrayOfObjects, HitCatalog &h1) {
	for (Object* obj: arrayOfObjects) {
		obj->rayTest(ray, h1);
	}
}

Vector3d refract_tray(Vector3d D, Vector3d N, double eta1, double eta2){
	double etar = eta1/eta2;
	double a = -1 *etar;
	double wn = D.dot(N);
	double radsq = pow(etar,2) * (pow(wn, 2) -1) +1;
	if (radsq < 0.0){
		Vector3d Ti(0,0,0);
		return Ti; 
	}
	else{
		double b = (etar * wn) - sqrt(radsq);
		Vector3d Ti = (a * D) + (b * N);
		return Ti;
	}
}
	
class Triangle: public Object{
	public:
		Vector3d point1;
		Vector3d point2;
		Vector3d point3;
		mtl materials;
		Vector3d N;
		Vector3d aMinusB;
		Vector3d aMinusC;
		double Ni;
		Vector3d vertexA;
		Vector3d vertexB;
		Vector3d vertexC;
		double sbeta; 
		double sgamma; 
		Triangle()=default;
		Triangle(Vector3d p1, Vector3d p2, Vector3d p3, mtl current){
			point1 = p1;
			point2 = p2;
			point3 = p3;
			materials = current;
			aMinusB = p1 - p2;
			aMinusC = p1 - p3;
		}
		Vector3d createN(HitCatalog h1){
			return N;
		}
		Ray refract_exit(Vector3d D, vector<Object*> arrayOfObjects, HitCatalog h1){
			Vector3d T1 = refract_tray(D, h1.N, 1, h1.mat.ni);
			Ray fakeVar; //this is completely a fake variable
			if(T1.sum() == 0.0){
				return fakeVar;
			}
			else{
				Vector3d exit; //this is a point
				Ray raye(h1.point, T1);
				HitCatalog h2;
				ray_find(raye, arrayOfObjects, h2);
				if(!h2.hit){
					return fakeVar;
				}
				Vector3d Nin = h2.N; 
				exit = h2.point;
				Nin.normalize();
				if (Nin.dot(raye.D) > 0){
						Nin *= -1;
				}
				Vector3d negT1 = -1*T1;
				Vector3d T2 = refract_tray(negT1, Nin, h1.mat.ni, 1);
				Ray refR(exit, T2);
				return refR;
			}
		}
		bool rayTest(Ray pixel, HitCatalog &h1) {
			Vector3d y = this->point1 - pixel.L;
			Matrix3d MM;
			MM(0, 0) = this->aMinusB[0];
			MM(1, 0) = this->aMinusB[1];
			MM(2, 0) = this->aMinusB[2];
			MM(0, 1) = this->aMinusC[0];
			MM(1, 1) = this->aMinusC[1];
			MM(2, 1) = this->aMinusC[2];
			MM(0, 2) = pixel.D[0];
			MM(1, 2) = pixel.D[1];
			MM(2, 2) = pixel.D[2];
			Matrix3d MMs1 = MM;
			Matrix3d MMs2 = MM;
			Matrix3d MMs3 = MM;
			for(int i=0;i<3;i++){
				MMs1(i,0) = y[i];
				MMs2(i,1) = y[i];
				MMs3(i,2) = y[i];
			}
			double detM  = MM.determinant();
			double detM1 = MMs1.determinant();
			double detM2 = MMs2.determinant();
			double sbeta = detM1 / detM; 
			double sgamma = detM2 / detM;
			if(sbeta>=0 and sgamma>=0 and sbeta+sgamma<=1){
				double detM3 = MMs3.determinant();
				double stval = detM3 / detM;
				if(stval>0.00001 and stval < h1.t){
					this->N = (1 - sbeta - sgamma) * this->vertexA + sbeta * this->vertexB + sgamma * this->vertexC;
					this->N.normalize();
					if (this->N.dot(pixel.D) > 0){
						this->N *= -1;
					}
					h1.N = this->N;
					h1.t = stval; 
					h1.point = pixel.L + (stval * pixel.D);
					h1.hit = true;
					h1.best = this;
					h1.mat = this->materials;
					return true;
				}
			}
			return false;
		}
};

class Sphere: public Object{
	public:
		Vector3d sphereCoords;
		Vector3d lightsReflected;
		double radius;
		Vector3d N;
		double Ni;
		Vector3d T;
		Vector3d Tr;
		mtl materials;
		Sphere()=default;
		Sphere(vector<string> data, int index){
			Vector3d coords(atof(data[index+1].c_str()), atof(data[index+2].c_str()), atof(data[index+3].c_str()));
			sphereCoords = coords;
			radius = atof(data[index+4].c_str());
			Vector3d ambient(atof(data[index+5].c_str()), atof(data[index+6].c_str()), atof(data[index+7].c_str()));
			Vector3d diffuse(atof(data[index+8].c_str()), atof(data[index+9].c_str()), atof(data[index+10].c_str()));
			Vector3d specular(atof(data[index+11].c_str()), atof(data[index+12].c_str()), atof(data[index+13].c_str()));
			Vector3d reflected(atof(data[index+14].c_str()), atof(data[index+15].c_str()), atof(data[index+16].c_str()));
			lightsReflected = reflected;
			Ni = atof(data[index+17].c_str());
			materials = mtl(16, specular, ambient, diffuse, reflected, Ni);
			Tr = materials.Tr;
		}
		Vector3d createN(HitCatalog h1){
			Vector3d N = h1.point - this->sphereCoords;
			N.normalize();
			return N;
		}
		Ray refract_exit(Vector3d D, vector<Object*> arrayOfObjects, HitCatalog h1){
			Ray fakeVar; //this variable is completely fake
			Vector3d T1 = refract_tray(D, h1.N, 1, h1.mat.ni);
			if(T1.sum() == 0.0){
				return fakeVar;
			}
			else{
				Vector3d exit = h1.point + 2*(this->sphereCoords-h1.point).dot(T1) *T1;
				Vector3d Nin = this->sphereCoords - exit;
				Nin.normalize();
				Vector3d negT1 = -1*T1;
				Vector3d T2 = refract_tray(negT1, Nin, h1.mat.ni, 1);
				Ray refR(exit, T2);
				return refR;
			}
		}
		bool rayTest(Ray pixel, HitCatalog &h1) { 
			Vector3d T = this->sphereCoords-pixel.L;
			double v    = T.dot(pixel.D);
			double csq  = T.dot(T);
			double disc = pow(this->radius,2) - (csq - pow(v,2)); 
			if (disc > 0.0) {
				double t  = v - sqrt(disc);
				if(t<h1.t and t>0.001){
					h1.t = t; //you might need to comment this out
					h1.point = pixel.L + (t * pixel.D);
					h1.hit = true;
					h1.best = this;
					h1.mat = this->materials;
					return true;
				}
			}
			return false;
		}
};

class Light{
	public:
		Vector3d coords;
		Vector3d rgb;
		double w;
	Light()=default;
	Light(vector<string> data,int &index){
		Vector3d lightCoords(atof(data[index+1].c_str()), atof(data[index+2].c_str()), atof(data[index+3].c_str()));
		w = atof(data[index+4].c_str());
		Vector3d lightRgb(atof(data[index+5].c_str()), atof(data[index+6].c_str()), atof(data[index+7].c_str()));
		coords = lightCoords;
		rgb = lightRgb;
	}
};

Ray pixelRay(int i, int j, Camera camera1){
	double px = ((i/(camera1.res["w"]-1.0))*(camera1.bounds["r"]-camera1.bounds["l"])+camera1.bounds["l"]);
	double py = ((j/(camera1.res["h"]-1.0))*(camera1.bounds["b"]-camera1.bounds["t"])+camera1.bounds["t"]);
	Vector3d L = camera1.eye + (camera1.nearFar["near"] * camera1.w) + (px * camera1.u) + (py * camera1.v);
	Vector3d D = L - camera1.eye;
	return Ray(L, D);
}

bool shadows(Vector3d pt, Light light, vector<Object*> arrayOfObjects){
	Vector3d toLight = light.coords - pt; //we are remaking this here because it was normalized earlier - we dont want that
	Ray ray(pt, toLight);
	double dtl = toLight.dot(ray.D);
	HitCatalog h1;
	ray_find(ray, arrayOfObjects, h1);
	if(h1.t < dtl){
		return true;
	}
	return false;
}

bool ray_trace(Ray ray, vector<Light> arrayOfLights, vector<Object*> arrayOfObjects, Vector3d ambient, Vector3d &accum, Vector3d refatt, int level) {
	Vector3d fake(-1000,-1000,-1000);
	HitCatalog h1;
	
	ray_find(ray, arrayOfObjects, h1);
	if (h1.hit){
		h1.N = h1.best->createN(h1); //N = surface normal
		Vector3d color = ambient.cwiseProduct(h1.mat.Ka);
		if (h1.N.dot(ray.D) > 0){
			h1.N *= -1;
		}
		for (int i=0; i<int(arrayOfLights.size());i++){
			Vector3d toLight = arrayOfLights[i].coords - h1.point;
			toLight.normalize();
			double NdotL = h1.N.dot(toLight);//lights specific stuff after this
			if (NdotL > 0.0 and not shadows(h1.point, arrayOfLights[i], arrayOfObjects)) {
				color += h1.mat.Kd.cwiseProduct(arrayOfLights[i].rgb) * NdotL; //to here - is lambertian component
				Vector3d toC = ray.L - h1.point; //ray back to camera - we are about to compute the far specular highlights
				toC.normalize();
				Vector3d spR = (2 * NdotL * h1.N) - toLight; //this is reflection
				spR.normalize();
				double CdR = toC.dot(spR);  //just dot product of ray back to camera and reflection
				if (CdR > 0.0) {   //if cosine of angle between them is postive, then you can take the material specular constants multiplied by the emission from the original light source, and multiply it by the cosine of hte angle between them raised to the alpha power (which is 16 for spheres)
					color += h1.mat.Ks.cwiseProduct(arrayOfLights[i].rgb)*pow(CdR, h1.mat.Ns);
				}
			}
		}
		for (int i=0;i<3;i++) { 
			accum[i] += refatt[i] * color[i]; // if all we were doing was recasting, we wouldve been done by now; we want to add color based on pretending that we are sitting on the surface looking at the actual reflection ray. means we need a new place for accumulating powers.
		}
		if (level > 0) { 
			Vector3d flec(0,0,0);
			Vector3d Uinv= -1* ray.D;
			Vector3d refR = (2 * h1.N.dot(Uinv) * h1.N) - Uinv; //refR
			refR.normalize();
			ray_trace(Ray(h1.point, refR), arrayOfLights, arrayOfObjects, ambient,flec, h1.mat.Kr.cwiseProduct(refatt), level-1);
			accum += flec;
		}
		if(level >0 and h1.mat.ni !=0 and h1.mat.Tr.sum() > 0){
			Vector3d thru(0,0,0);
			Vector3d negD = -1 * ray.D;
			Ray fraR = h1.best->refract_exit(negD, arrayOfObjects, h1);
			if (fraR.D != fake){
				ray_trace(fraR, arrayOfLights, arrayOfObjects, ambient, thru, h1.mat.Tr.cwiseProduct(refatt), level-1);
				accum += thru;
			} 
		} 
	}
	return true;
}
	
int main(int argc, char *argv[]){ //read in driver file
	int frames = 1;
	string s(argv[1]);
	if(s== "animation"){
		system("mkdir animationppms");
		frames = 180;
	}
	for(int fi=0;fi<frames;fi++){
		ifstream read;
		ofstream write;
		string drName; 
		if(argc ==2){
			drName = "driver";
			if(to_string(fi).length()<2){
				drName += "00" + to_string(fi);
			}
			else if(to_string(fi).length()<3){
				drName += "0" + to_string(fi);
			}
			else{
				drName += to_string(fi);
			}
			read.open(drName + ".txt");
		}
		else if(argc == 3){
			read.open(argv[1]);
		}
		else {
			return 0;
		}
		string token;
		vector<string> data;
		string driver = argv[1];
		while(!read.eof()){
			read >> token;
			if(!read.eof()){
			data.push_back(token);
			}
			else{
				break;
			}
		}
		vector<Object*> arrayOfObjects;
		vector<Light> arrayOfLights;
		int recursionLevel = 0;
		double focalpoint = numeric_limits<double>::max();
		double apsize = 1;
		int rays = 1;
		Vector3d ambient;
		Camera camera1;
		for(int i=0;i<int(data.size());i++){ //camera
			if(data[i]=="focalpoint"){
				focalpoint = atof(data[i+1].c_str());
				apsize = atof(data[i+2].c_str());
				rays = atof(data[i+3].c_str());
			}
			if(data[i]=="recursionlevel"){
				recursionLevel = atof(data[i+1].c_str());
			}
			if(data[i]=="eye"){
				Vector3d e(atof(data[i+1].c_str()), atof(data[i+2].c_str()), atof(data[i+3].c_str()));
				camera1.eye = e;
			}
			if(data[i]=="look"){
				Vector3d l(atof(data[i+1].c_str()), atof(data[i+2].c_str()), atof(data[i+3].c_str()));
				camera1.look = l;
			}	
			if(data[i]=="up"){
				Vector3d upp(atof(data[i+1].c_str()), atof(data[i+2].c_str()), atof(data[i+3].c_str()));
				camera1.up = upp;
			}
			if(data[i]=="d"){
				camera1.nearFar["near"] = atof(data[i+1].c_str());
			}
			if(data[i]=="bounds"){
				camera1.bounds["l"] = atof(data[i+1].c_str());
				camera1.bounds["r"] = atof(data[i+2].c_str());
				camera1.bounds["b"] = atof(data[i+3].c_str());
				camera1.bounds["t"] = atof(data[i+4].c_str());
			}
			if(data[i]=="res"){
				camera1.res["w"] = atof(data[i+1].c_str());
				camera1.res["h"] = atof(data[i+2].c_str());
			}
			if(data[i]=="ambient"){
				Vector3d amb(atof(data[i+1].c_str()), atof(data[i+2].c_str()), atof(data[i+3].c_str())); 
				ambient = amb;
			}
			if(data[i]=="light"){
				Light ligh = Light(data, i);
				arrayOfLights.push_back(ligh);
			}
			if(data[i]=="sphere"){
				Sphere* sph = new Sphere(data, i);
				arrayOfObjects.push_back(sph);
			}
			if(data[i] == "model"){
				vector<string> models;
				
				//rotate
				//rotate step 1 - normalize vector w
				Vector3d w(atof(data[i+1].c_str()),atof(data[i+2].c_str()),atof(data[i+3].c_str()));
				w.normalize();
				int minimum = 0;
				//rotate step 2 - find min value, set equal to 1, create new vector m
				if(w[0] <= w[1] && w[0] <= w[2]){
					minimum = 0;
				}
				else if(w[1] <= w[2] && w[1] <= w[0]){
					minimum = 1;
				}
				else if(w[2] <= w[1] && w[2] <= w[0]){
					minimum = 2;
				}
				Vector3d m = w;
				m[minimum] = 1;
				m.normalize();
				//rotate step 3 - create vector u
				Vector3d u = w.cross(m);
				u.normalize();
				//rotate step 4 - create vector v
				Vector3d v = w.cross(u);
				
				Matrix4d R;
				R << 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					0,0,0,1;
				R(0,0) = u[0];
				R(0,1) = u[1];
				R(0,2) = u[2];
				R(1,0) = v[0];
				R(1,1) = v[1];
				R(1,2) = v[2];
				R(2,0) = w[0];
				R(2,1) = w[1];
				R(2,2) = w[2];
				
				//read in object file
				ifstream object(data[i+10]);
				string word;
				mtl current;
				while(word!="mtllib"){
					object >> word;
				}
				string mtlname;
				object >> mtlname;
				
				ifstream mtlfile(mtlname);
				string tok;
				mtlfile >> tok;
				map<string,mtl> mtlArray;
				while(tok!="newmtl"){
					mtlfile >> tok;
				}
				while(!mtlfile.eof()){
					if(tok=="newmtl"){
						mtl newmtl;
						mtlfile >> tok;
						newmtl.name = tok;
						bool TrExist = false;
						while(tok!="newmtl" and !mtlfile.eof()){
							mtlfile >> tok;
							if(tok == "Ns"){
								mtlfile >> tok;
								newmtl.Ns = atof(tok.c_str());
							}
							if(tok == "Ka"){
								string tok1;
								string tok2;
								string tok3;
								mtlfile >> tok1;
								mtlfile >> tok2;
								mtlfile >> tok3;
								Vector3d K(atof(tok1.c_str()), atof(tok2.c_str()), atof(tok3.c_str()));
								newmtl.Ka = K;
							}
							if(tok == "Kd"){
								string tok1;
								string tok2;
								string tok3;
								mtlfile >> tok1;
								mtlfile >> tok2;
								mtlfile >> tok3;
								Vector3d K(atof(tok1.c_str()), atof(tok2.c_str()), atof(tok3.c_str()));
								newmtl.Kd = K;
							}
							if(tok == "Ks"){
								string tok1;
								string tok2;
								string tok3;
								mtlfile >> tok1;
								mtlfile >> tok2;
								mtlfile >> tok3;
								Vector3d K(atof(tok1.c_str()), atof(tok2.c_str()), atof(tok3.c_str()));
								newmtl.Ks = K;
							}
							if(tok == "Ke"){
								string tok1;
								string tok2;
								string tok3;
								mtlfile >> tok1;
								mtlfile >> tok2;
								mtlfile >> tok3;
								Vector3d K(atof(tok1.c_str()), atof(tok2.c_str()), atof(tok3.c_str()));
								newmtl.Ke = K;
							}
							if(tok == "Ni"){
								mtlfile >> tok;
								newmtl.ni = atof(tok.c_str());
							}
							if(tok == "Tr"){
								string tok1;
								string tok2;
								string tok3;
								mtlfile >> tok1;
								mtlfile >> tok2;
								mtlfile >> tok3;
								Vector3d K(atof(tok1.c_str()), atof(tok2.c_str()), atof(tok3.c_str()));
								newmtl.Tr = K;
								TrExist = true;
							}
							if(tok == "d"){
								mtlfile >> tok;
								newmtl.di = atof(tok.c_str());
							}
							if(tok == "illum"){
								mtlfile >> tok;
								newmtl.illum = atof(tok.c_str());
							}
						}
						Vector3d zeroes(0,0,0);
						if(!TrExist){
							newmtl.Tr = zeroes;
						}
						if(newmtl.illum ==2){
							newmtl.Kr = zeroes;
						}
						else{
							newmtl.Kr = newmtl.Ks;
						}
						mtlArray[newmtl.name]=newmtl;
						
					}
				}
				object >> word;
				int q = 0;
				map<int, mtl> mtlMap;
				if(word=="usemtl"){
					object >> word;
					current = mtlArray[word];
					object >> word;
					mtlMap[q] = current;
					q++;
				}
				int j = 0;
				MatrixXd obj(4,j);
				vector<vector<int>> faces;
				map<int, vector<vector<int>>> vNes;
				while(!object.eof()){
					if(word == "v"){
						obj.conservativeResize(4,j+1);
						object >> word;
						obj(0,j) = atof(word.c_str());
						object >> word;
						obj(1,j) = atof(word.c_str());
						object >> word;
						obj(2,j) = atof(word.c_str());
						obj(3,j) = 1;
						j++;
					}
				//Making triangles
					if(word=="usemtl"){
						object >> word;
						current = mtlArray[word];
						mtlMap[q] = current;
						q++;
					}
					if(word=="f"){
						vector<int> fs;
						object >> word;
						fs.push_back(atoi(word.substr(0,word.length()-3).c_str())-1);
						object >> word;
						fs.push_back(atoi(word.substr(0,word.length()-3).c_str())-1);
						object >> word;
						fs.push_back(atoi(word.substr(0,word.length()-3).c_str())-1);
						fs.push_back(q-1);
						vNes[fs[0]].push_back(fs);
						vNes[fs[1]].push_back(fs);
						vNes[fs[2]].push_back(fs);
						faces.push_back(fs);
					}
					object >> word;
				}
				MatrixXd transformation(4, obj.cols());
				
				//create rotation matrix 
				Matrix4d rotate;
				rotate << 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					0,0,0,1;
				double angle = atof(data[i+4].c_str());
				angle = angle*M_PI/180;
				rotate(0,0) = cos(angle);
				rotate(1,1) = cos(angle);
				rotate(0,1) = -sin(angle);
				rotate(1,0) = sin(angle);
				Matrix4d Rt = R.inverse();
				//complete rotation, and do scaling at same time
				Matrix4d scale;
				scale << 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					0,0,0,1;
				scale(0,0) =atof(data[i+5].c_str());
				scale(1,1) =atof(data[i+5].c_str());
				scale(2,2) = atof(data[i+5].c_str());
				transformation = scale*Rt*rotate*R;
				Matrix4d p;
				p << 1,0,0,0,
					0,1,0,0,
					0,0,1,0,
					0,0,0,1;
				p(0,3) = atof(data[i+6].c_str());
				p(1,3) = atof(data[i+7].c_str());
				p(2,3) = atof(data[i+8].c_str());
				
				transformation = p*transformation;
				
				MatrixXd transformedObj(4, obj.cols());
				transformedObj = transformation*obj;
				
				//Vertices
				vector<Vector3d> vertices;
				for(int j=0; j<transformedObj.cols();j++){
					Vector3d vertice(transformedObj(0,j), transformedObj(1,j), transformedObj(2,j));
					vertices.push_back(vertice);
				}
				double cutoff = atof(data[i+9].c_str());
				for(int p =0; p<int(faces.size());p++){
					vector<int> fs = faces[p];
					Triangle* tri = new Triangle(vertices[fs[0]], vertices[fs[1]], vertices[fs[2]], mtlMap[fs[3]]);
					Vector3d line1 = tri->point2 - tri->point1;
					Vector3d line2 = tri->point3 - tri->point1;
					Vector3d N = line1.cross(line2);
					N.normalize();
					Vector3d none(0,0,0);
					tri->vertexA = none;
					tri->vertexB = none;
					tri->vertexC = none;
					for(int whichVertex = 0; whichVertex < 3; whichVertex++){
						vector<vector<int>> faceList = vNes[fs[whichVertex]];
						for (int f2 = 0;f2<int(faceList.size());f2++){
							line1 = vertices[faceList[f2][1]] - vertices[faceList[f2][0]];
							line2 = vertices[faceList[f2][2]] - vertices[faceList[f2][0]];
							Vector3d tempN = line1.cross(line2);
							tempN.normalize();
							double NdotN = N.dot(tempN);
							if (NdotN>1){
								NdotN = 1;
							}
							if (NdotN<-1){
								NdotN = -1;
							}
							if(int(acos(NdotN)*180/M_PI)<=cutoff){
								switch(whichVertex){
									case 0: tri->vertexA = tri->vertexA + tempN;
										break;
									case 1: tri->vertexB = tri->vertexB + tempN;
										break;
									case 2: tri->vertexC = tri->vertexC + tempN;
										break;
								}
							}
						}
						
					}
					tri->vertexA = tri->vertexA / tri->vertexA.norm();
					tri->vertexB = tri->vertexB / tri->vertexB.norm();
					tri->vertexC = tri->vertexC / tri->vertexC.norm();
					arrayOfObjects.push_back(tri);
				} 
			}
		}   
		camera1.makeCoords();
		ofstream newfile;
		string newName;
		if(argc ==2){
			newName =  drName+ ".ppm";
		}
		else{
			newName = argv[2];
		}
		newfile.open(newName);
		newfile << "P3\n"<< camera1.res["w"] << " " << camera1.res["h"] << " 255\n";
		//Ray pixel;
		auto p =[](double x) -> double { return max(0.0,min(255.0,floor(255.0 * x)));};
		vector<vector<Vector3d>> allColors (camera1.res["w"]);
		for (int i =0;i< camera1.res["w"];i++){
			allColors[i].resize(camera1.res["h"]);
		}
		int width = camera1.res["w"];
		int height = camera1.res["h"];
		#pragma omp parallel for collapse(2)
		for (int i=0; i<width;i++){
			for (int j=0; j< height;j++){
				Ray pixel = pixelRay(j, i, camera1);
				Vector3d totColor(0.0,0.0,0.0);
				if(rays == 1){
					Vector3d refatt(1.0,1.0,1.0);
					int level = recursionLevel;
					ray_trace(pixel, arrayOfLights, arrayOfObjects, ambient, totColor, refatt, level);
				}
				else{
					Vector3d P = camera1.eye + focalpoint * pixel.D;
					Vector3d firstRefatt(1.0,1.0,1.0);
					Vector3d firstColor(0.0,0.0,0.0);
					int firstLevel = recursionLevel;
					ray_trace(pixel, arrayOfLights, arrayOfObjects, ambient, firstColor, firstRefatt, firstLevel);
					totColor+=firstColor;
					for (int k = 0;k<rays;k++){
						Vector3d color(0.0,0.0,0.0);
						double rw = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * 2.0 - 1.0;
						double rh = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * 2.0 - 1.0;
						double rz = (static_cast<float>(rand() % RAND_MAX) / RAND_MAX) * 2.0 - 1.0;
						Vector3d aperture((apsize*rw), (apsize*rh), (apsize*rz));
						Vector3d apertureSubtract(abs(camera1.w[0])*aperture[0], abs(camera1.w[1])*aperture[1], abs(camera1.w[2])*aperture[2]);
						aperture = aperture - apertureSubtract;
						Vector3d L = camera1.eye + aperture;
						Ray ray(L, P-L);
						Vector3d refatt(1.0,1.0,1.0);
						int level = recursionLevel;
						ray_trace(ray, arrayOfLights, arrayOfObjects, ambient, color, refatt, level);
						totColor += color;
					}
					totColor = totColor/(rays+1);
				}
				allColors[i][j] = totColor;
			}
		}
		for(int w=0;w<camera1.res["w"];w++){
			for (int h=0;h<camera1.res["h"];h++){
				newfile << p(allColors[w][h][0]) << " " << p(allColors[w][h][1]) << " " << p(allColors[w][h][2]) << " ";
			}
			newfile << "\n";
		}
		newfile.close();
		if(argc ==2){
			system(("mv " + newName + " animationppms").c_str());
		}
	}
	return 0;
}
