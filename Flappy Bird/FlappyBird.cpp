#include "GlobalVariable.h"

#define DEG2RAD 3.1459/180.0
float horizontal, vertical;
float tickf;
float yLocation = 0.0f; // Keep track of our position on the y axis.  
bool movingUp = false; // Whether or not we are moving up or down  
float yRotationAngle = 0.0f; // The angle of rotation for our object  
int orange_egg = 0, white_egg = 0, pink_egg = 0, yellow_egg = 0, gray_egg = 0, green_egg = 0;
int score;
char sScore[100];
char test[50];
int egg_xc, egg_yc;
int d = 0; //for color selection
float rotation;

char flagX,flagY; 
float tick[8],xFunc[8]={-700,-700-200,-700-400,-700-600,-700-800,-700-1000,-700-1200},yFunc,speed=0.7;
float YPipeRoof,YPipeButtom,Ybird,threshYdown[8],threshYup[8];

typedef struct {
	float m[4][4];
} matrix3D_t;

typedef struct {
	float v[4];
} vector3D_t;

typedef struct {
	float x;
	float y;
	float z;
} point3D_t;

typedef struct {
	float x;
	float y;
} point2D_t;

typedef struct {
	float r;
	float g;
	float b;
} color_t;

//////////////////SETTING/////////////////////////////////////

////////////////// matrices and vectors 3D ver 2 /////////////////
matrix3D_t createIdentity(void)
{
	matrix3D_t u;
	int i, j;
	for (i = 0; i<4; i++) {
		for (j = 0; j<4; j++) u.m[i][j] = 0.;
		u.m[i][i] = 1.;
	}
	return u;
}

matrix3D_t operator * (matrix3D_t a, matrix3D_t b)
{
	matrix3D_t c;//c=a*b
	int i, j, k;
	for (i = 0; i<4; i++) for (j = 0; j<4; j++) {
		c.m[i][j] = 0;
		for (k = 0; k<4; k++) c.m[i][j] += a.m[i][k] * b.m[k][j];
	}
	return c;
}

vector3D_t operator * (matrix3D_t a, vector3D_t b)
{
	vector3D_t c;//c=a*b
	int i, j;
	for (i = 0; i<4; i++) {
		c.v[i] = 0;
		for (j = 0; j<4; j++) c.v[i] += a.m[i][j] * b.v[j];
	}
	return c;
}

matrix3D_t translationMTX(float dx, float dy, float dz)
{
	matrix3D_t trans = createIdentity();
	trans.m[0][3] = dx;
	trans.m[1][3] = dy;
	trans.m[2][3] = dz;
	return trans;
}

matrix3D_t rotationXMTX(float theta)
{
	matrix3D_t rotate = createIdentity();
	float cs = cos(theta);
	float sn = sin(theta);
	rotate.m[1][1] = cs; rotate.m[1][2] = -sn;
	rotate.m[2][1] = sn; rotate.m[2][2] = cs;
	return rotate;
}

matrix3D_t rotationYMTX(float theta)
{
	matrix3D_t rotate = createIdentity();
	float cs = cos(theta);
	float sn = sin(theta);
	rotate.m[0][0] = cs; rotate.m[0][2] = sn;
	rotate.m[2][0] = -sn; rotate.m[2][2] = cs;
	return rotate;
}

matrix3D_t rotationZMTX(float theta)
{
	matrix3D_t rotate = createIdentity();
	float cs = cos(theta);
	float sn = sin(theta);
	rotate.m[0][0] = cs; rotate.m[0][1] = -sn;
	rotate.m[1][0] = sn; rotate.m[1][1] = cs;
	return rotate;
}

matrix3D_t scalingMTX(float factorx, float factory, float factorz)
{
	matrix3D_t scale = createIdentity();
	scale.m[0][0] = factorx;
	scale.m[1][1] = factory;
	scale.m[2][2] = factorz;
	return scale;
}

matrix3D_t perspectiveMTX(float eyelength)
{
	matrix3D_t perspective = createIdentity();
	perspective.m[3][2] = -1. / eyelength;
	return perspective;
}

point2D_t Vector2Point2D(vector3D_t vec)
{
	point2D_t pnt;
	pnt.x = vec.v[0];
	pnt.y = vec.v[1];
	return pnt;
}

point3D_t Vector2Point3D(vector3D_t vec)
{
	point3D_t pnt;
	pnt.x = vec.v[0];
	pnt.y = vec.v[1];
	pnt.z = vec.v[2];
	return pnt;
}

vector3D_t Point2Vector(point3D_t pnt)
{
	vector3D_t vec;
	vec.v[0] = pnt.x;
	vec.v[1] = pnt.y;
	vec.v[2] = pnt.z;
	vec.v[3] = 1.;
	return vec;
}

vector3D_t homogenizeVector(vector3D_t vec)
{
	int i;
	for (i = 0; i<3; i++) {
		vec.v[i] /= vec.v[3];
	}
	vec.v[3] = 1.;
	return vec;
}

vector3D_t unitVector(vector3D_t vec)
{
	int i;
	float vec2 = 0.;
	float vec1, invvec1;
	for (i = 0; i<3; i++) {
		vec2 += vec.v[i] * vec.v[i];
	}
	vec1 = sqrt(vec2);
	if (vec1 != 0.) {
		invvec1 = 1. / vec1;
		for (i = 0; i<3; i++) {
			vec.v[i] *= invvec1;
		}
	}
	vec.v[3] = 1.;
	return vec;
}

// inner product (dot product) of homogeneous vector
float operator * (vector3D_t a, vector3D_t b)
{
	float c;//c=a*b
	int i;
	c = 0;
	for (i = 0; i<3; i++) {
		c += a.v[i] * b.v[i];
	}
	return c;
}

// outer product (cross product ) of homogeneous vector
//       i         j         k
//       a0       a1        a2
//       b0       b1        b2
vector3D_t operator ^ (vector3D_t a, vector3D_t b)
{
	vector3D_t c;//c=a*b
	c.v[0] = a.v[1] * b.v[2] - a.v[2] * b.v[1];
	c.v[1] = a.v[2] * b.v[0] - a.v[0] * b.v[2];
	c.v[2] = a.v[0] * b.v[1] - a.v[1] * b.v[0];
	c.v[3] = 1.;
	return c;
}

vector3D_t operator - (vector3D_t v1, vector3D_t v0)
{
	vector3D_t c;//c=v1-v0
	c.v[0] = v1.v[0] - v0.v[0];
	c.v[1] = v1.v[1] - v0.v[1];
	c.v[2] = v1.v[2] - v0.v[2];
	c.v[3] = 1.;
	return c;
}

vector3D_t operator - (vector3D_t v)
{
	vector3D_t c;//c=-v
	c.v[0] = -v.v[0];
	c.v[1] = -v.v[1];
	c.v[2] = -v.v[2];
	c.v[3] = 1.;
	return c;
}

vector3D_t operator * (float r, vector3D_t b)
{
	vector3D_t c;//c=r*b
	int i;
	for (i = 0; i<3; i++) {
		c.v[i] = r*b.v[i];
	}
	c.v[3] = 1.;
	return c;
}

vector3D_t operator * (vector3D_t b, float r)
{
	vector3D_t c;//c=r*b
	int i;
	for (i = 0; i<3; i++) {
		c.v[i] = r*b.v[i];
	}
	c.v[3] = 1.;
	return c;
}

float funcPositive(float x)
{
	if (0.<x) return x;
	else return 0.;
}

// x to yth power
float power(float x, float y)
{
	//ln z = y ln x        z = exp (y ln x)
	if (x == 0.) return 0;
	return exp(y*log(x));
}

color_t operator + (color_t c1, color_t c2)
{
	color_t col;
	col.r = c1.r + c2.r;
	col.g = c1.g + c2.g;
	col.b = c1.b + c2.b;
	return col;
}

color_t operator * (float r, color_t c)
{
	color_t col;
	col.r = r*c.r;
	col.g = r*c.g;
	col.b = r*c.b;
	return col;
}

color_t operator * (color_t c, float r)
{
	color_t col;
	col.r = r*c.r;
	col.g = r*c.g;
	col.b = r*c.b;
	return col;
}

//PhongModel color calculation
// LightVector, NormalVector, ViewVector, ColorofObject
color_t PhongModel(vector3D_t Light, vector3D_t Normal, vector3D_t View, color_t col)
{
	float kspe = 0.7; // specular reflection coefficient
	float kdif = 0.6; // diffuse reflection coefficient
	float kamb = 0.4; // ambient light coefficient
	float tmp, NL, RV;
	color_t ColWhite = { 1, 1, 1 };
	vector3D_t ReflectionVector = (2.*(Light*Normal)*Normal) - Light;
	tmp = Normal*Light;
	NL = funcPositive(tmp);
	tmp = ReflectionVector*View;
	RV = funcPositive(tmp);
	return kdif*NL*col + kspe*power(RV, 4)*ColWhite + kamb*col;
}

////////////// End of matrices and vectors 3D ver 2 //////////////
////////////// OpenGL drawShape Functions ver 1 /////////////////
void setColor(float red, float green, float blue)
{
	glColor3f(red, green, blue);
}

void setColor(color_t col)
{
	glColor3f(col.r, col.g, col.b);
}

void drawDot(point2D_t p)
{
	glBegin(GL_POINTS);
	glVertex2f(p.x, p.y);
	glEnd();
}

void drawLine(float x1, float y1, float x2, float y2)
{
	glBegin(GL_LINES);
	glVertex2f(x1, y1);
	glVertex2f(x2, y2);
	glEnd();
}

void drawLine(point2D_t p1, point2D_t p2)
{
	drawLine(p1.x, p1.y, p2.x, p2.y);
}

//n: number of points
void drawPolyline(point2D_t pnt[], int n)
{
	int i;
	glBegin(GL_LINE_STRIP);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

//n: number of vertices
void drawPolygon(point2D_t pnt[], int n)
{
	int i;
	glBegin(GL_POLYGON);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

// The function fillPolygon can fills only convex polygons
//n: number of vertices
void fillPolygon(point2D_t pnt[], int n, color_t color)
{
	int i;
	setColor(color);
	glBegin(GL_POLYGON);
	for (i = 0; i<n; i++) {
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

// The function gradatePolygon can fills only convex polygons
// The vertices will be painted with corresponding given colors.
// The points inside the polygon will be painted with the mixed color.
//n: number of vertices
void gradatePolygon(point2D_t pnt[], int num, color_t col[])
{
	int i;
	glBegin(GL_POLYGON);
	for (i = 0; i<num; i++) {
		setColor(col[i]);
		glVertex2f(pnt[i].x, pnt[i].y);
	}
	glEnd();
}

//////////// End of OpenGL drawShape Functions ver 1 ////////////

////////////// OpenGL drawShape Functions ver 1 /////////////////
typedef struct {
	int NumberofVertices; //in the face
	short int pnt[50];
	color_t col;
} face_t;
typedef struct {
	int NumberofVertices; //of the object
	point3D_t pnt[1600];
	color_t col[1600];
	int NumberofFaces; //of the object
	face_t fc[1000];
} object3D_t;

void colors()
{
	switch (d)
	{
	case 3:glColor3f(1, 0.4, 0); break; //orange
	case 2:glColor3f(0, 1, 0); break; //green
	case 4:glColor3f(0.2, 0.2, 0.2); break; //gray
	case 5:glColor3f(1, 1, 1); break; //white
	case 1:glColor3f(1, 1, 0); break; //yellow
	}
}

void draw3D(object3D_t obyek, matrix3D_t mat){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	point2D_t p[50];
	int i, j;
	for (i = 0; i<obyek.NumberofVertices; i++){
		vec[i] = Point2Vector(obyek.pnt[i]);
		vec[i] = mat*vec[i];
	}
	setColor(1, 1, 0);
	for (i = 0; i<obyek.NumberofFaces; i++){
		for (j = 0; j<obyek.fc[i].NumberofVertices; j++)
			vecbuff[j] = vec[obyek.fc[i].pnt[j]];
		vecNormal = (vecbuff[1] - vecbuff[0]) ^ (vecbuff[2] - vecbuff[0]);
		if (vecNormal.v[2]<0){
			for (j = 0; j<obyek.fc[i].NumberofVertices; j++){
				p[j] = Vector2Point2D(vecbuff[j]);
			}
			drawPolygon(p, obyek.fc[i].NumberofVertices);
		}

	}
	setColor(1, 1, 0);
	colors();
	for (i = 0; i<obyek.NumberofFaces; i++){
		for (j = 0; j<obyek.fc[i].NumberofVertices; j++)
			vecbuff[j] = vec[obyek.fc[i].pnt[j]];
		vecNormal = (vecbuff[1] - vecbuff[0]) ^ (vecbuff[2] - vecbuff[0]);
		if (vecNormal.v[2] >= 0){
			for (j = 0; j<obyek.fc[i].NumberofVertices; j++){
				p[j] = Vector2Point2D(vecbuff[j]);
			}
			drawPolygon(p, obyek.fc[i].NumberofVertices);

		}

	}

}

void createCube(object3D_t &kubus, float d){
	object3D_t obyek = { 8,
	{ { 0, 0, 0 }, { d, 0, 0 }, { d, 0, d }, { 0, 0, d }, { 0, d, 0 }, { d, d, 0 }, { d, d, d }, { 0, d, d } },
	{ { 1, 1, 0 }, { 1, 1, 0 }, { 1, 1, 0 }, { 1, 1, 0 }, { 1, 1, 0 }, { 1, 1, 0 }, { 1, 1, 0 } },
	6,
	{ { 4, { 0, 1, 2, 3 }, { 1, 1, 0 } }, { 4, { 4, 7, 6, 5 }, { 1, 1, 0 } }, { 4, { 1, 5, 6, 2 }, { 1, 1, 0 } }, { 4, { 0, 3, 7, 4 }, { 1, 1, 0 } }, { 4, { 2, 6, 7, 3 }, { 1, 1, 0 } }, { 4, { 0, 4, 5, 1 }, { 1, 1, 0 } } } };
	matrix3D_t mat = translationMTX(-d / 2, -d / 2, -d / 2);
	vector3D_t vec;
	kubus = obyek;
	for (int i = 0; i<kubus.NumberofVertices; i++){
		vec = Point2Vector(kubus.pnt[i]);
		vec = mat*vec;
		kubus.pnt[i] = Vector2Point3D(vec);
	}
}
void drawDot3D(point3D_t p, matrix3D_t mat, color_t col){
	vector3D_t vec, lightVector = { 0, 0, 1, 1 }, viewVector = { 0, 0, 1, 1 };
	point2D_t pt;
	color_t colbuff;
	vec = Point2Vector(p);
	vec = mat*vec;
	pt = Vector2Point2D(vec);
	glPointSize(4);
	vec = unitVector(vec);
	colbuff = PhongModel(lightVector, vec, viewVector, col);
	setColor(colbuff);
	drawDot(pt);
}
void createSphere(object3D_t &sphere, int n, float r){
	float a = 6.28 / n;
	float b = 6.28 / n;
	int i, j;
	sphere.NumberofVertices = (n + 1)*n;
	for (i = 0; i <= n; i++){
		for (j = 0; j<n; j++){
			sphere.pnt[i*n + j].x = r*cos(j*a)*sin(i*b);
			sphere.pnt[i*n + j].y = r*cos(i*b);
			sphere.pnt[i*n + j].z = r*sin(j*a)*sin(i*b);
		}
	}
	sphere.NumberofFaces = n*n + 2;
	for (i = 0; i<n; i++){
		for (j = 0; j<n; j++){
			sphere.fc[i*n + j].NumberofVertices = 4;
			sphere.fc[i*n + j].pnt[0] = i*n + j;
			sphere.fc[i*n + j].pnt[1] = (i + 1)*n + j;
			sphere.fc[i*n + j].pnt[2] = (i + 1)*n + j + 1;
			sphere.fc[i*n + j].pnt[3] = i*n + j + 1;
			if (j == (n - 1)){
				sphere.fc[i*n + j].pnt[2] = i*n + j + 1;
				sphere.fc[i*n + j].pnt[3] = (i - 1)*n + j + 1;
			}
		}
	}
	sphere.fc[n*n].NumberofVertices = n;
	for (i = 0; i<n; i++) sphere.fc[n*n].pnt[i] = i;
	sphere.fc[n*n + 1].NumberofVertices = n;
	for (i = 0; i<n; i++)
		sphere.fc[n*n + 1].pnt[i] = (n + 1)*n - 1 - i;
	color_t c = { 1, 1, 0 };
	for (i = 0; i<sphere.NumberofFaces; i++)
		sphere.fc[i].col = c;
	for (i = 0; i<sphere.NumberofVertices; i++)
		sphere.col[i] = c;
}

void createCylinderN(object3D_t &silinder, int m, int n, float r[], float
	h[]){
	float a = 6.26 / n;
	float b = 0;
	int i, j;
	silinder.NumberofVertices = (m + 1)*n;
	for (i = 0; i <= m; i++){
		if (i>0) b = b + h[i - 1];
		for (j = 0; j<n; j++){
			silinder.pnt[i*n + j].x = r[i] * cos(j*a);
			silinder.pnt[i*n + j].y = b;
			silinder.pnt[i*n + j].z = r[i] * sin(j*a);
		}
	}
	silinder.NumberofFaces = m*n + 2;
	for (i = 0; i<m; i++){
		for (j = 0; j<n; j++){
			silinder.fc[i*n + j].NumberofVertices = 4;
			silinder.fc[i*n + j].pnt[0] = i*n + j;
			silinder.fc[i*n + j].pnt[1] = (i + 1)*n + j;
			silinder.fc[i*n + j].pnt[2] = (i + 1)*n + j + 1;
			silinder.fc[i*n + j].pnt[3] = i*n + j + 1;
			if (j == (n - 1)){
				silinder.fc[i*n + j].pnt[2] = i*n + j + 1;
				silinder.fc[i*n + j].pnt[3] = (i - 1)*n + j + 1;
			}
		}
	}
	silinder.fc[m*n].NumberofVertices = n;
	for (i = 0; i<n; i++) silinder.fc[m*n].pnt[i] = i;
	silinder.fc[m*n + 1].NumberofVertices = n;
	for (i = 0; i<n; i++)
		silinder.fc[m*n + 1].pnt[i] = (m + 1)*n - 1 - i;
}

//============================code NEW


color_t DiffuseScattering(vector3D_t Light, vector3D_t Normal,
	vector3D_t View, color_t col) {
	// diffuse reflection coefficient 
	float kdif = 0.6;
	float tmp, NL, RV;
	color_t ColWhite = { 1, 1, 1 };
	vector3D_t ReflectionVector = (2.*(Light*Normal)*Normal) - Light;
	tmp = Normal*Light;
	NL = funcPositive(tmp);
	tmp = ReflectionVector*View;
	RV = funcPositive(tmp);
	return kdif*NL*col;
}

color_t SpecuralReflection(vector3D_t Light, vector3D_t Normal,
	vector3D_t View, color_t col) {
	// specular reflection coefficient 
	float kspe = 0.7;
	float tmp, NL, RV;
	color_t ColWhite = { 1, 1, 1 };
	vector3D_t ReflectionVector = (2.*(Light*Normal)*Normal) - Light;
	tmp = Normal*Light;
	NL = funcPositive(tmp);
	tmp = ReflectionVector*View;
	RV = funcPositive(tmp);
	return kspe*power(RV, 4)*ColWhite;
}

color_t Ambient(vector3D_t Light, vector3D_t Normal, vector3D_t View,
	color_t col) {
	// ambient light coefficient 
	float kamb = 0.4;
	return kamb*col;
}

color_t LambertLaw(vector3D_t Light, vector3D_t Normal, vector3D_t
	View, color_t col) {
	// diffuse reflection coefficient 
	float kdif = 0.6;
	// ambient light coefficient 
	float kamb = 0.4;
	float tmp, NL, RV;
	color_t ColWhite = { 1, 1, 1 };
	vector3D_t ReflectionVector = (2.*(Light*Normal)*Normal) - Light;
	tmp = Normal*Light;
	NL = funcPositive(tmp);
	tmp = ReflectionVector*View;
	RV = funcPositive(tmp);
	return kdif*NL*col + kamb*col;
}




void draw3Dc(object3D_t obyek, matrix3D_t mat, color_t col){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	vector3D_t lightVector = { 0, 0, 1, 1 };
	vector3D_t viewVector = { 0, 0, 1, 1 };
	color_t colbuff;
	point2D_t p[50];
	int i, j;
	for (i = 0; i < obyek.NumberofVertices; i++){
		vec[i] = Point2Vector(obyek.pnt[i]);
		vec[i] = mat*vec[i];
	}
	for (i = 0; i < obyek.NumberofFaces; i++){
		for (j = 0; j < obyek.fc[i].NumberofVertices; j++)
			vecbuff[j] = vec[obyek.fc[i].pnt[j]];
		vecNormal = (vecbuff[1] - vecbuff[0]) ^ (vecbuff[2] - vecbuff[0]);
		if (vecNormal.v[2] < 0){
			for (j = 0; j < obyek.fc[i].NumberofVertices; j++){
				p[j] = Vector2Point2D(vecbuff[j]);
			}
			vecNormal = unitVector(vecNormal);
			colbuff = PhongModel(lightVector, vecNormal, viewVector, col);
			fillPolygon(p, obyek.fc[i].NumberofVertices, colbuff);
		}

	}
	for (i = 0; i < obyek.NumberofFaces; i++){
		for (j = 0; j < obyek.fc[i].NumberofVertices; j++)
			vecbuff[j] = vec[obyek.fc[i].pnt[j]];
		vecNormal = (vecbuff[1] - vecbuff[0]) ^ (vecbuff[2] - vecbuff[0]);
		if (vecNormal.v[2] >= 0){
			for (j = 0; j < obyek.fc[i].NumberofVertices; j++){
				p[j] = Vector2Point2D(vecbuff[j]);
			}
			vecNormal = unitVector(vecNormal);
			colbuff = PhongModel(lightVector, vecNormal, viewVector, col);
			fillPolygon(p, obyek.fc[i].NumberofVertices, colbuff);
		}
	}
}

color_t warna={0,1,0};

void Pipa(float x,float tinggi,float lebar,char direction){
	float pembagi;float y;
	if(tinggi<100){pembagi= 3;}
	else pembagi=4;

	

if(direction==1){
	y=-0.24*x+(-156);//Funtion of Y
	glLoadIdentity();
	glTranslatef(x, y, 0.0f);//atur posisi
	matrix3D_t tilting = rotationXMTX(0.25)*rotationYMTX(-0.5 + yLocation);
	object3D_t tabung;
	
	float r[4] = { lebar, lebar, lebar, lebar };//atur lebar
	float h[3] = { tinggi, tinggi, tinggi };//atur tinggi
	createCylinderN(tabung, 3, 20, r, h);
	draw3Dc(tabung, tilting,warna);
    
	glLoadIdentity();
	YPipeButtom = y+tinggi*2.35;
	glTranslatef(x,YPipeButtom, 0.0f);//atur posisi
	float r2[4] = { lebar+lebar*0.3, lebar+lebar*0.3, lebar+lebar*0.3, lebar+lebar*0.3 };//atur lebar
	float h2[3] = { tinggi/pembagi, tinggi/pembagi, tinggi/pembagi};//atur tinggi
	createCylinderN(tabung, 3, 20, r2, h2);
	draw3Dc(tabung, tilting,warna);
}

else if(direction==0){

	y=-0.24*x+(-156)+700;//Funtion of Y
	glLoadIdentity();
	glTranslatef(x, y, 0.0f);//atur posisi
	matrix3D_t tilting = rotationXMTX(0.25)*rotationYMTX(-0.5 + yLocation);
	object3D_t tabung;
	
	float r[4] = { lebar, lebar, lebar, lebar };//atur lebar
	float h[3] = { tinggi, tinggi, tinggi };//atur tinggi
	createCylinderN(tabung, 3, 20, r, h);
	draw3Dc(tabung, tilting,warna);
    
	glLoadIdentity();
	YPipeRoof   = y+tinggi*2.8;
	glTranslatef(x,YPipeRoof, 0.0f);//atur posisi
	float r2[4] = { lebar+lebar*0.3, lebar+lebar*0.3, lebar+lebar*0.3, lebar+lebar*0.3 };//atur lebar
	float h2[3] = { tinggi/10, tinggi/10, tinggi/10};//atur tinggi
	createCylinderN(tabung, 3, 20, r2, h2);
	draw3Dc(tabung, tilting,warna);

}





}



void gridLines(){


	for (int i = -320; i <= 320; i = i + 20){

		glLineWidth(2);
		glBegin(GL_LINES);
		if (i == 0){ glColor3f(1.0, 0.0, 0.0); }
		else glColor3f(1.0, 1.0, 1.0);
		glVertex2f(-650, i);
		glVertex2f(650, i);
		glEnd();
		glFlush();

	}
	for (int i = -320; i <= 320; i = i + 20){


		glLineWidth(2);
		glBegin(GL_LINES);
		if (i == 0){ glColor3f(1.0, 0.0, 0.0); }
		else glColor3f(1.0, 1.0, 1.0);
		glVertex2f(i, -350);
		glVertex2f(i, 350);
		glEnd();
		glFlush();

	}

}

void drawEllipse(float cx, float cy, float radiusX, float radiusY){

	glBegin(GL_POLYGON);glColor3f(1.0, 0.0, 1.0);
	for (int i = 0; i < 360; i++){

		float rad = i*DEG2RAD;

		glVertex2f(cos(rad)*radiusX + cx, sin(rad)*radiusY + cy);
	}
	glEnd();
	glFlush();

}


void Write(char text[], int x, int y, float r, float g, float b){
	glColor4f(r, g, b, 0.0f);

	glRasterPos2f(x, y);

	const char * p = text;
	do glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *p);
	while (*(++p));



}

//==========================================================================================================//
//==========================================================================================================//
//==========================================================================================================//
//==========================================================================================================//

//=================	menu=======================//

static int window, returnmenu, usemenu, value = 0;

void menu(int n){
	if (n == 0)
	{
		glutDestroyWindow(window);
		exit(0);
	}

	else
	{
		value = n;
	}
	glutPostRedisplay();
}
void createMenu(void){
	// membuat sub menu
	usemenu = glutCreateMenu(menu);
	//glutAddMenuEntry("Play Flappy Bird", 2);

	//membuat menu
	returnmenu = glutCreateMenu(menu);
	//glutAddSubMenu("Play!", usemenu);
	glutAddMenuEntry("Play!",2);
	glutAddMenuEntry("Quit", 0);
	
	//mengaktifkan fungsi button kanan dari mouse
	glutAttachMenu(GLUT_RIGHT_BUTTON);

}
//===========================================//






void background(){
		// langit
		glBegin(GL_POLYGON);
		glColor3f(0.556862, 0.737254, 0.925490);
		glVertex3f(-650, -350, 0);
		glVertex3f(-650, 350, 0);
		glVertex3f(650, 350, 0);
		glVertex3f(650, -350, 0);
		glEnd();
		glFlush();

		// tanah
		glBegin(GL_POLYGON);
		glColor3f(0.509803, 0.368627, 0.282352);
		glVertex3f(-650, -350, 0);
		glVertex3f(-650, 50, 0);
		glVertex3f(650, -200, 0);
		glVertex3f(650, -350, 0);
		glEnd();
		glFlush();

		// awan
		glColor3f(1, 1, 1);
		glBegin(GL_POLYGON);
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + 360, sin(rad) * 20 + 275);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 25 + 400, sin(rad) * 22 + 275);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + 440, sin(rad) * 20 + 275);
		}
		glEnd();
		glFlush();

		// awan
		glColor3f(1, 1, 1);
		glBegin(GL_POLYGON);
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + 410, sin(rad) * 20 + 225);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 25 + 450, sin(rad) * 22 + 225);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + 490, sin(rad) * 20 + 225);
		}
		glEnd();
		glFlush();

		// awan
		glColor3f(1, 1, 1);
		glBegin(GL_POLYGON);
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + -360, sin(rad) * 20 + 275);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 25 + -400, sin(rad) * 22 + 275);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + -440, sin(rad) * 20 + 275);
		}
		glEnd();
		glFlush();

		// awan
		glColor3f(1, 1, 1);
		glBegin(GL_POLYGON);
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + -410, sin(rad) * 20 + 225);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 25 + -450, sin(rad) * 22 + 225);
		}
		for (int i = 0; i < 360; i++){
			float rad = i * DEG2RAD;
			glVertex2f(cos(rad) * 30 + -490, sin(rad) * 20 + 225);
		}
		glEnd();
		glFlush();
	}

//=======================================FOR BIRD==================================//

void body(){
	//matrix3D_t tilting = rotationXMTX(0.25)*rotationYMTX(-0.5);
	//matrix3D_t tilting = rotationXMTX(0.25 + rotation)*rotationYMTX(-0.5);
	matrix3D_t tilting = rotationXMTX(0)*rotationYMTX(-0.25 );
	setColor(0, 1, 0);
	object3D_t bola;
	createSphere(bola, 20, 30);
	setColor(0, 0, 0);
	color_t col = { 1, 1, 0 };
	draw3Dc(bola, tilting,col);
}

void mata(){
	//matrix3D_t tilting = rotationXMTX(0.25)*rotationYMTX(-0.5);
	//matrix3D_t tilting = rotationXMTX(0.25 + rotation)*rotationYMTX(-0.5);
	matrix3D_t tilting = rotationXMTX(0)*rotationYMTX(-0.25);
	setColor(0, 1, 0);
	object3D_t bola;
	createSphere(bola, 20, 10);
	setColor(0, 0, 0);
	color_t col = { 1, 1, 1 };
	draw3Dc(bola, tilting, col); 
}

void dalammata(){
	//matrix3D_t tilting = rotationXMTX(0.25)*rotationYMTX(-0.5);
	//matrix3D_t tilting = rotationXMTX(0.25 + rotation)*rotationYMTX(-0.5);
	matrix3D_t tilting = rotationXMTX(0)*rotationYMTX(-0.25 );
	setColor(0, 1, 0);
	object3D_t bola;
	createSphere(bola, 20, 5);
	setColor(0, 0, 0);
	color_t col = { 0, 0, 0 };
	draw3Dc(bola, tilting, col);
}

void paruh(){
	rotation -= 9;
	//matrix3D_t tilting = rotationXMTX(0)*rotationYMTX(-0.25 );
	matrix3D_t tilting = rotationXMTX(0)*rotationYMTX(-0.25 );
	setColor(0, 1, 0);
	object3D_t o;
	float r[10] = { 0, 0, 0, 25, 30, 30, 25, 0, 0, 0 };
	float h[19] = { 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0 };
	createCylinderN(o, 9, 20, r, h);
	color_t col = { 1, 0, 0 };
	draw3Dc(o, tilting, col);

}

void paruh1(){
	rotation -= 9;
	matrix3D_t tilting = rotationXMTX(0 )*rotationYMTX(-0.25);
	setColor(0, 1, 0);
	object3D_t o;
	float r[10] = { 0, 0, 0, 25, 30, 30, 25, 0, 0, 0 };
	float h[19] = { 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2, 1, 0, 0 };
	createCylinderN(o, 9, 20, r, h);
	color_t col = { 1, 0, 0 };
	draw3Dc(o, tilting,col);
 

}
  
float xbird = 50;

void Bird(){

	
	
    static int dir,L =0,flag=0;static float tick=-200;

	if(flagJump==1){
		Ybird = tickf+=5.5;
	    drawEllipse(xbird,Ybird,15,15);
		if(tickf>=180){tickf=180;}
    }
    else {Ybird = tickf-=2;
    	drawEllipse(xbird,Ybird,15,15);
		if(tickf<=-200){tickf=-200;}}

}


void BirdNew(){ readJoystick();
static int tickBird;
	if(flagJump==1){Ybird= tickBird+=5; if(tickBird>250)tickBird=250;}
	else if(flagJump==0){Ybird = tickBird-=2; if(tickBird<=-180)tickBird=-180;}

	glPushMatrix();
	glTranslatef(-10, -15+tickBird, 0);
	paruh();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-10, -20+tickBird, 0);
	paruh1();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-25, 10+tickBird, 0);
	mata();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(-27.5, 10+tickBird, 0);
	dalammata();
	glPopMatrix();
	glPushMatrix();
	glTranslatef(0, 0+tickBird, 0);
	body();
	glPopMatrix();
}






//================================END BIRD==================================//
void TimerInterrupt(int timer){
readJoystick();
 glutTimerFunc(1/1000,TimerInterrupt,1);//1ms
}
void pipaRUN(){

	Pipa(xFunc[0],20,15+tick[0],1);threshYdown[0]=YPipeButtom;Pipa(xFunc[0],-150,15+tick[0],0);threshYup[0]=YPipeRoof;//Pipa 1
    Pipa(xFunc[1],40,15+tick[1],1);threshYdown[1]=YPipeButtom;Pipa(xFunc[1],-120,15+tick[1],0);threshYup[1]=YPipeRoof;//Pipa 2
    Pipa(xFunc[2],30,15+tick[2],1);threshYdown[2]=YPipeButtom;Pipa(xFunc[2],-150,15+tick[2],0);threshYup[2]=YPipeRoof;//Pipa 3
    Pipa(xFunc[3],60,15+tick[3],1);threshYdown[3]=YPipeButtom;Pipa(xFunc[3],-100,15+tick[3],0);threshYup[3]=YPipeRoof;//Pipa 4
    Pipa(xFunc[4],30,15+tick[4],1);threshYdown[4]=YPipeButtom;Pipa(xFunc[4],-170,15+tick[4],0);threshYup[4]=YPipeRoof;//Pipa 5
    Pipa(xFunc[5],40,15+tick[5],1);threshYdown[5]=YPipeButtom;Pipa(xFunc[5],-150,15+tick[5],0);threshYup[5]=YPipeRoof;//Pipa 6
    Pipa(xFunc[6],50,15+tick[6],1);threshYdown[6]=YPipeButtom;Pipa(xFunc[6],-120,15+tick[6],0);threshYup[6]=YPipeRoof;//Pipa 7
  
     //readJoystick();
	if(xFunc[0]>700){xFunc[0]=-700;}	else {xFunc[0]+=speed;}
	if(xFunc[1]>700){xFunc[1]=-700;}  	else {xFunc[1]+=speed;}
	if(xFunc[2]>700){xFunc[2]=-700;}	else {xFunc[2]+=speed;}
	if(xFunc[3]>700){xFunc[3]=-700;}	else {xFunc[3]+=speed;}
	if(xFunc[4]>700){xFunc[4]=-700;}	else {xFunc[4]+=speed;}
	if(xFunc[5]>700){xFunc[5]=-700;}	else {xFunc[5]+=speed;}
	if(xFunc[6]>700){xFunc[6]=-700;}	else {xFunc[6]+=speed;}
	// readJoystick();

	if(xFunc[0]>-650){tick[0]+=0.012;}	else {tick[0]=0;}
	if(xFunc[1]>-650){tick[1]+=0.012;}  else {tick[1]=0;}
	if(xFunc[2]>-650){tick[2]+=0.012;}	else {tick[2]=0;}
	if(xFunc[3]>-650){tick[3]+=0.012;}	else {tick[3]=0;}
	if(xFunc[4]>-650){tick[4]+=0.012;}	else {tick[4]=0;}
	if(xFunc[5]>-650){tick[5]+=0.012;}	else {tick[5]=0;}
	if(xFunc[6]>-650){tick[6]+=0.012;}	else {tick[6]=0;}
 	

	//printf("flag.jump = %d\n",flagJump);



}

void Aiflappy(void){
	
static int step=1;int yPipeBawah,yPipeAtas,xPipa;
	if(step==1){yPipeBawah=threshYdown[0];yPipeAtas=threshYup[0];xPipa=xFunc[0];}
	else if(step==2){yPipeBawah=threshYdown[1];yPipeAtas=threshYup[1];xPipa=xFunc[1];}
	else if(step==3){yPipeBawah=threshYdown[2];yPipeAtas=threshYup[2];xPipa=xFunc[2];}
	else if(step==4){yPipeBawah=threshYdown[3];yPipeAtas=threshYup[3];xPipa=xFunc[3];}
	else if(step==5){yPipeBawah=threshYdown[4];yPipeAtas=threshYup[4];xPipa=xFunc[4];}
	else if(step==6){yPipeBawah=threshYdown[5];yPipeAtas=threshYup[5];xPipa=xFunc[5];}
	else if(step==7){yPipeBawah=threshYdown[6];yPipeAtas=threshYup[6];xPipa=xFunc[6];}


	if(Ybird<yPipeBawah || Ybird>yPipeAtas){
		//printf("masuk kondisi Y\n");
		if(xbird==xPipa+60){value=3;}
	
	}
	else {
		if(xbird==xPipa){step++;if(step>7)step=1;}

	}
	printf("%d\n",step);

}

void userdraw(void){
	
	
	   //gridLines();
		glBegin(GL_LINES);
		/*glColor3f(1.0, 0.0, 0.0);
		glVertex2f(-650, 0);
		glVertex2f(650, 0);
		glColor3f(1.0, 0.0, 0.0); 
		glVertex2f(0, -350);
		glVertex2f(0, 350);
		*/
		/*glColor3f(0.0, 1.0, 1.0);
		glVertex2f(-650, 50);
		glVertex2f(650, 50);
*/
		////glColor3f(0.0, 0.5, 1.0);
		////glVertex2f(-650, 0);
		////glVertex2f(600, -300);

/*
		glColor3f(0.0, 0.5, 1.0);
		glVertex2f(-650, 0+220);
		glVertex2f(600, -300+220);
*/
/*
		glColor3f(0.0, 0.5, 1.0);
		glVertex2f(-650, 100);
		glVertex2f(600, 400);
*/
		glEnd();
		glFlush();


		Write("CALIBRATE YOUR CONSOLE",0,0,1,1,1);
	    glFlush();

		
		if(value==1)return ; 
		else if (value==2){
			Aiflappy();
			background();
			//Bird();
 			BirdNew();
			pipaRUN();
		}
		else if (value == 3) {Write("GameOver",0,0,1,1,1);}
	

	//tick[0]+=0.004;if(tick[0]>10){tick[0]=10;}

}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glLoadIdentity();
	userdraw();
	glutSwapBuffers();
}

void reshape(int w, int h) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-650., 650., -350.0, 350.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0, 0, w, h);
}

int main(int argc, char**argv)
{
	//=========Init Serial===============
	signal(SIGINT, signalHandler);
	signal(SIGTERM, signalTerminated);
	system(SERIAL_SYS);
	serial.serialPort = openSerialPort();
	sleep(1);
	//===================================




	glutInit (&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(1300, 700);
	glutCreateWindow("Teknik Komputer");
	glClearColor(0.0, 0.0, 0.0, 0.0);
	gluOrtho2D(-650., 650., -350.0, 350.0);
	createMenu();
	glutIdleFunc(display);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape); // Tell GLUT to use the method "reshape" for reshaping  
	glutTimerFunc(1/1000,TimerInterrupt,1);
	glutMainLoop();
	return 0;
}