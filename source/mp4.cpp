
#include <GLEW/glew.h>
#include <SOIL.h>
#include <GLUT/glut.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

#define PI 3.1415926535

//wireframe on vs off
bool wireframe = false;

// the x,y,z coordinates of the I
double xco[] = {-1,-.18,.18,1,-1,-.18,.18,1,-1,-.18,.18,1,-1,-.18,.18,1,-1,-.18,.18,1,-1,-.18,.18,1,-1,-.18,.2,1,-1,-.18,.18,1};

double yco[] = {1.2,1.2,1.2,1.2,.75,.75,.75,.75,-.75,-.75,-.75,-.75,-1.2,-1.2,-1.2,-1.2,1.2,1.2,1.2,1.2,.75,.75,.75,.75,-.75,-.75,-.75,-.75,-1.2,-1.2,-1.2,-1.2};

double zco[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5,-.5};

//vertices that make the faces
int facelist[] = {
5, 6, 2, 1, 4, 5, 1, 0, 6, 7, 3, 2, 9, 10, 6, 5, 12, 13, 9, 8, 13, 14, 10, 9, 14, 15, 
11, 10, 16, 17, 21, 20,17, 18, 22, 21, 18, 19, 23, 22, 21, 22, 26, 25, 24, 25, 29, 28, 25, 26, 30, 29, 26,
 27, 31, 30, 17, 16, 0, 1, 18, 17, 1, 2, 19, 18, 2, 3, 13, 12, 28, 29, 14, 13, 29, 30, 15, 14, 30, 31, 5, 4, 
 20, 21, 7, 6, 22, 23, 8, 9, 25, 24, 10, 11, 27, 26, 3, 7, 23, 19,11, 15, 31, 27,6, 10, 26, 22, 16, 20, 4, 0,21, 25, 9, 5, 24, 28, 12, 8}; 


//To hold the mesh we keep track of all vertices, all faces and all of the half edges

class mesh_hedge;

//the vertex class to hold a vertex
class mesh_vertex{
public:
	//default constructor
	mesh_vertex(){
		vertex_he = NULL;
		x = 0; y = 0; z = 0;
		nx = 0; ny = 0; nz = 0;
	}
	//constructor takes in 3 doubles
	mesh_vertex(double xc,double yc,double zc){
		vertex_he = NULL;
		x = xc; y = yc; z = zc;
		nx = 0; ny = 0; nz = 0;
	}
	//a half egde of the vertex
	mesh_hedge * vertex_he;
	//the position of the vertex
	double x, y, z;
	//the vertex normal
	double nx, ny, nz;
};
//the face class to hold a face
class mesh_face{
public:
	//constructors for the mesh
	mesh_face(){
		face_he = NULL;
		face_vertex = NULL;
		a = b = c = d = NULL;
	}
	mesh_face(mesh_vertex * ax,mesh_vertex * bx,mesh_vertex * cx,mesh_vertex * dx){
		face_he = NULL;
		face_vertex = NULL;
		a = ax; b = bx; c = cx; d = dx;
	}
	mesh_face(mesh_vertex * ax,mesh_vertex * bx,mesh_vertex * cx,mesh_vertex * dx,mesh_hedge * e){
		face_he = e;
		face_vertex = NULL;
		a = ax; b = bx; c = cx; d = dx;
	}
	//the four vertices of the face
	mesh_vertex *a;
	mesh_vertex *b;
	mesh_vertex *c; 
	mesh_vertex *d;
	//the face vertex
	mesh_vertex * face_vertex;
	//a halfegde connected to the face
	mesh_hedge * face_he; 
	
};
//half edge class
class mesh_hedge{
public:
	//constructors
	mesh_hedge(){
		vertex_a = NULL;vertex_b = NULL;vertex_e = NULL;
		he_h0 = NULL;he_h1 = NULL;he_n = NULL;he_p = NULL;
		face_he = NULL;
	}
	mesh_hedge(mesh_vertex * s,mesh_vertex * e,mesh_face * f){
		he_n = NULL;he_p = NULL;vertex_e = NULL;he_h1 = NULL;he_h0 = NULL;
		vertex_a = s;
		vertex_b = e;
		face_he = f;
	}
	//variables of the half edge
	mesh_vertex * vertex_a;
	mesh_vertex * vertex_b;
	mesh_vertex* vertex_e; 
	mesh_hedge* he_h0;
	mesh_hedge* he_h1;
	mesh_hedge * he_n; 
	mesh_hedge * he_p;
	mesh_face * face_he; 
};


//vectors to hold the mesh
vector<mesh_hedge*> edges;
vector<mesh_vertex*> vertices;
vector<mesh_face*> faces;

//calculate the face vertices by going through each face and averaging all 4 vertices 
void calcfaces(){
	for(int i = 0; i < faces.size(); i++){
		faces[i]->face_vertex = new mesh_vertex((faces[i]->a->x + faces[i]->b->x + faces[i]->c->x + faces[i]->d->x)/4,(faces[i]->a->y + faces[i]->b->y + faces[i]->c->y + faces[i]->d->y)/4,(faces[i]->a->z + faces[i]->b->z + faces[i]->c->z + faces[i]->d->z)/4);
		vertices.push_back(faces[i]->face_vertex);
	}
}
//calculate the edge vertices by looping through the edges
void calcedges(){

	for(int i =0;i<edges.size();i++){
		if( edges[i]->he_p->vertex_e != NULL)
			edges[i]->vertex_e = edges[i]->he_p->vertex_e;
		else{
			edges[i]->vertex_e = new mesh_vertex((edges[i]->face_he->face_vertex->x + edges[i]->he_p->face_he->face_vertex->x + edges[i]->vertex_b->x + edges[i]->he_p->vertex_b->x)/4,(edges[i]->face_he->face_vertex->y + edges[i]->he_p->face_he->face_vertex->y + edges[i]->vertex_b->y + edges[i]->he_p->vertex_b->y)/4,(edges[i]->face_he->face_vertex->z + edges[i]->he_p->face_he->face_vertex->z + edges[i]->vertex_b->z + edges[i]->he_p->vertex_b->z)/4);
			vertices.push_back(edges[i]->vertex_e);
		}
	}
}
vector<mesh_vertex*> verts;
//calculate the new vertices by looping through them
void calcverts(){

	vector<mesh_hedge*> hearray;
	vector<mesh_vertex*> vfaces;
	vector<mesh_vertex*> vedges;
	
	for(int i=0; i < verts.size(); i++)
	{
		hearray.clear();
		vfaces.clear();
		vedges.clear();

		hearray.push_back(verts[i]->vertex_he);
		for(int r=1;r<5;r++)
			hearray.push_back(hearray[r-1]->he_n->he_n->he_n->he_p);

		for(int r=0; r<5;r++){
			vfaces.push_back(hearray[r]->vertex_e);
			vedges.push_back(hearray[r]->face_he->face_vertex);
		}

		double num,m,num2;

		if(hearray[0] == hearray[3]){
			num=3.0;m=1.0;num2=4.0;	
		}
		else if(hearray[0] == hearray[4]){
			num=4.0;m=1.0;num2=4.0;	
		}
		else{
			num=5.0;m=2.0;num2=5.0;	
		}

		vfaces.push_back(new mesh_vertex());
		vedges.push_back(new mesh_vertex());

		//take the average
		for(int r=0;r<num;r++){
			vfaces[5]->x+=vfaces[r]->x/num;
			vfaces[5]->y+=vfaces[r]->y/num;
			vfaces[5]->z+=vfaces[r]->z/num;

			vedges[5]->x+=vedges[r]->x/num;
			vedges[5]->y+=vedges[r]->y/num;
			vedges[5]->z+=vedges[r]->z/num;
		}
		//set the new vertex corordinates
		verts[i]->x = (m*verts[i]->x +  2*vedges[5]->x + vfaces[5]->x)/num2; 
		verts[i]->y = (m*verts[i]->y +  2*vedges[5]->y + vfaces[5]->y)/num2; 
		verts[i]->z = (m*verts[i]->z +  2*vedges[5]->z + vfaces[5]->z)/num2; 	
	}
}

vector<mesh_face*> tempfaces;
vector<mesh_hedge*> tempedges;

//creates the mesh using all of new vertices
void createmesh(){
	//clear out the temp vars
	tempfaces.clear();
	tempedges.clear();
	vector<mesh_hedge*> hearray;
	vector<mesh_face*> facearray;
		
	for(int i = 0; i<faces.size(); i++){
		
		hearray.clear();
		facearray.clear();
		for(int r=0;r<20;r++){
			hearray.push_back(new mesh_hedge());

		}	

		hearray[0]  = faces[i]->face_he;
		for(int r=1;r<4;r++)
			hearray[r]  = hearray[r-1]->he_n;
		//set up all the half edges
		for(int r=0;r<4;r++){

			hearray[r]->he_h1 =hearray[4+2*r];
			hearray[r]->he_h0 = hearray[5+2*r];	
			hearray[4+2*r]->vertex_a = hearray[r]->vertex_a;
			hearray[4+2*r]->vertex_b = hearray[r]->vertex_e;
			hearray[5+2*r]->vertex_a = hearray[r]->vertex_e;
			hearray[5+2*r]->vertex_b = hearray[r]->vertex_b; 
			if(hearray[r]->he_p->he_h1 != NULL){
				hearray[5+2*r]->he_p = hearray[r]->he_p->he_h1;
				hearray[4+2*r]->he_p = hearray[r]->he_p->he_h0;
				
				hearray[r]->he_p->he_h1->he_p = hearray[5+2*r];
				hearray[r]->he_p->he_h0->he_p = hearray[4+2*r];
			}

			hearray[r]->vertex_a->vertex_he = hearray[4+2*r];

			hearray[12+2*r]->vertex_a = hearray[r]->vertex_e;
			hearray[12+2*r]->vertex_b = faces[i]->face_vertex;
			hearray[13+2*r]->vertex_b = hearray[r]->vertex_e;
			hearray[13+2*r]->vertex_a = faces[i]->face_vertex;
			hearray[12+2*r]->he_p = hearray[13+2*r];
			hearray[13+2*r]->he_p = hearray[12+2*r];

			hearray[r]->vertex_e->vertex_he = hearray[12+2*r];

			hearray[4+2*r]->he_n = hearray[12+2*r];
		}

		faces[i]->face_vertex->vertex_he = hearray[19];
	
		hearray[12]->he_n = hearray[19];
		hearray[19]->he_n = hearray[11]; 
		hearray[11]->he_n = hearray[4];
		
		for(int r=0;r<3;r++){
			hearray[14+2*r]->he_n = hearray[13+2*r];
			hearray[13+2*r]->he_n = hearray[5+2*r]; 
			hearray[5+2*r]->he_n = hearray[6+2*r];
		}
		//create the new faces	
		facearray.push_back(new mesh_face(hearray[4]->vertex_a,hearray[12]->vertex_a,hearray[19]->vertex_a,hearray[11]->vertex_a,hearray[11]));
		facearray.push_back(new mesh_face(hearray[6]->vertex_a,hearray[14]->vertex_a,hearray[13]->vertex_a,hearray[5]->vertex_a,hearray[5]));
		facearray.push_back(new mesh_face(hearray[8]->vertex_a,hearray[16]->vertex_a,hearray[15]->vertex_a,hearray[7]->vertex_a,hearray[7]));
		facearray.push_back(new mesh_face(hearray[10]->vertex_a,hearray[18]->vertex_a,hearray[17]->vertex_a,hearray[9]->vertex_a,hearray[9]));
		
		int indices[]= {0,1,1,2,2,3,3,0,0,1,1,2,2,3,3,0};
		//set the half edge faces
		for(int r = 0; r<16;r++)
		{
			hearray[4+r]->face_he = facearray[indices[r]];
		}
		//add the faces and egdes
		for(int r=0;r<4;r++)
			tempfaces.push_back(facearray[r]);
		for(int r = 4; r<20;r++)
			tempedges.push_back(hearray[r]);
	}
}

//struct to hold a normal vector
struct nvec{
	double x,y,z;
};

//this methed goes though every mesh_vertex in the mesh and computes the per mesh_vertex normals
void calcnormals(){

	vector<mesh_hedge*> hearray;
	nvec avgNorm;
	vector<nvec> normals;

	
	for(int i = 0; i < vertices.size(); i ++){
		
		hearray.clear();
		//get all the half edges
		hearray.push_back(vertices[i]->vertex_he);
		for(int r=1;r<5;r++)
			hearray.push_back(hearray[r-1]->he_n->he_n->he_n->he_p);

		normals.clear();
		//get all the normals
		for(int r=0;r<5;r++){
			nvec normal;
			
			mesh_vertex a = *(hearray[r]->vertex_a);
			mesh_vertex b = *(hearray[r]->vertex_b);
			mesh_vertex c = *(hearray[r]->he_n->vertex_b);

			normal.x = ((b.y-a.y)*(c.z-b.z)) - ((c.y-b.y)*(b.z-a.z));
			normal.y = ((b.z-a.z)*(c.x-b.x)) - ((c.z-b.z)*(b.x-a.x));
			normal.z = ((b.x-a.x)*(c.y-b.y)) - ((c.x-b.x)*(b.y-a.y));

			normals.push_back(normal);
		}
		nvec t;
		normals.push_back(t);

		double num;
		if(hearray[0] == hearray[3])
			num=3.0;
		else if(hearray[0] == hearray[4])
			num=4.0;
		else
			num=5.0;
		//average the normals
		for(int r=0;r<num;r++){
			normals[5].x+=normals[r].x/num;
			normals[5].y+=normals[r].y/num;
			normals[5].z+=normals[r].z/num;

		}
		//make sure length of normal vector is one
		float length = sqrt(normals[5].x*normals[5].x + normals[5].y*normals[5].y + normals[5].z*normals[5].z);

		normals[5].x /=length;
		normals[5].y /=length;
		normals[5].z /=length;
		//set the normal
		vertices[i]->nx = normals[5].x;
		vertices[i]->ny = normals[5].y;
		vertices[i]->nz = normals[5].z;

	}
}
//clears out all the mesh structures and then addes the unsubdived I to the mesh
void clearmesh(){
	
	//clear the  vectors holding the mesh
	faces.clear();
	edges.clear();
	vertices.clear();
	//add all the initial vertices
	for(int i =0;i<32;i++){
		mesh_vertex * v = new mesh_vertex;
		v->x =xco[i];
		v->y =yco[i];
		v->z =zco[i];
		vertices.push_back(v);
	}
	//set up all the faces
	for(int i=0;i<120;i+=4){
	//cout<<"The i val is:"<<i<<vertex_bl;
	//cout<<"The mesh_face val is:"<<facelist[i]<<endlfacelist;

		mesh_face *cmesh_face = new mesh_face(vertices[facelist[i]],vertices[facelist[i+1]],vertices[facelist[i+2]],vertices[facelist[i+3]]);

		mesh_hedge *e1 = new mesh_hedge(vertices[facelist[i]],vertices[facelist[i+1]],cmesh_face);
		mesh_hedge *e2 = new mesh_hedge(vertices[facelist[i+1]],vertices[facelist[i+2]],cmesh_face);
		mesh_hedge *e3 = new mesh_hedge(vertices[facelist[i+2]],vertices[facelist[i+3]],cmesh_face);
		mesh_hedge *e4 = new mesh_hedge(vertices[facelist[i+3]],vertices[facelist[i]],cmesh_face);
		vertices[facelist[i]]->vertex_he = e1;
		vertices[facelist[i+1]]->vertex_he = e2;
		vertices[facelist[i+2]]->vertex_he = e3;
		vertices[facelist[i+3]]->vertex_he = e4;

		e1->he_n = e2;
		e2->he_n = e3;
		e3->he_n = e4;
		e4->he_n = e1;

		cmesh_face->face_he = e4;

		faces.push_back(cmesh_face);
		edges.push_back(e1);
		edges.push_back(e2);
		edges.push_back(e3);
		edges.push_back(e4);	

	}

	for( int i = 0; i < edges.size(); i++)
	{
		for(int j = i +1; j < edges.size(); j++)
		{
			if(edges[i]->vertex_a == edges[j]->vertex_b && edges[i]->vertex_b == edges[j]->vertex_a){
				edges[i]->he_p = edges[j];
				edges[j]->he_p = edges[i];
			}

		}
	}


}
//vars for the bezier curve
int direction;
double camtime;
//initialize the program
void init(void) 
{	
	direction = 1;
	camtime = 1.0;
	glewInit();

	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	
	//set the clear color
	glClearColor (1.0, 1.0, 1.0, 1.0);

	glEnable(GL_TEXTURE_2D);
	
	//load the textures
	GLuint texture;
	glGenTextures(1, &texture);
	
	int img_width, img_height;
	
	//Code to texture the I, Taken from mp3

	//use soil to load the image
	unsigned char* textimg = SOIL_load_image("stones.jpg", &img_width, &img_height, NULL, 0);
	
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, img_width, img_height, 0, GL_RGB, GL_UNSIGNED_BYTE, textimg);

	//Enable line smoothing so that the wireframe looks better
	glEnable( GL_LINE_SMOOTH );
	glEnable( GL_POLYGON_SMOOTH );
	glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
	glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );

	glLineWidth(.003);

	//initialize lighing
	GLfloat spec[] = {1.0,1.0,1.0};
	
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT0, GL_AMBIENT, spec);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, spec);
	glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

	//set up mesh to the intital I
	clearmesh();
	calcnormals();

}

void display(void)
{
	glLoadIdentity ();   
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//repeat bezier curve
	camtime += direction*.003;
	if(camtime >=1 || camtime<=0)
		direction = - direction;

	//look at the center of the I, camera moves on a bezier curve
	gluLookAt(
		pow((1-camtime),3)*1.78 + 3*pow((1-camtime),2)*camtime*-.67 + 3*(1-camtime)*pow(camtime,2)*.46 +  pow(camtime,3)*.46, 
		pow((1-camtime),3)*-1.78 + 3*pow((1-camtime),2)*camtime*3.89+ 3*(1-camtime)*pow(camtime,2)*5.89 +  pow(camtime,3)*-3.4, 
		pow((1-camtime),3)*.15 + 3*pow((1-camtime),2)*camtime*4.52 + 3*(1-camtime)*pow(camtime,2)*-.46 +  pow(camtime,3)*-3.4, 
		0,-.25,0,0,1,0);

	//Draw the mesh by looping through each mesh_face
	for(int i=0;i<faces.size();i++){	

		//if we want to show just the wireframe use a line loop
		if(wireframe){
			glBegin(GL_LINE_LOOP);
				glVertex3f(faces[i]->a->x,faces[i]->a->y, faces[i]->a->z); 
				glVertex3f(faces[i]->b->x,faces[i]->b->y, faces[i]->b->z); 
				glVertex3f(faces[i]->c->x,faces[i]->c->y, faces[i]->c->z); 
				glVertex3f(faces[i]->d->x,faces[i]->d->y, faces[i]->d->z); 
			glEnd();
		}
		//otherwise we use GL_QUADS
		else{
			glBegin(GL_QUADS);
			//calculate the texture coordinates using a cylinder
			glNormal3f(faces[i]->a->nx,faces[i]->a->ny,faces[i]->a->nz);
			glTexCoord2f((atan2(faces[i]->a->z +.25, faces[i]->a->x) + PI)/(2*PI), (faces[i]->a->y + 1.2)/2.4);
			glVertex3f(faces[i]->a->x,faces[i]->a->y, faces[i]->a->z); 

			glNormal3f(faces[i]->b->nx,faces[i]->b->ny,faces[i]->b->nz);
			glTexCoord2f((atan2(faces[i]->b->z +.25, faces[i]->b->x) + PI)/(2*PI), (faces[i]->b->y + 1.2)/2.4);
			glVertex3f(faces[i]->b->x,faces[i]->b->y, faces[i]->b->z); 

			glNormal3f(faces[i]->c->nx,faces[i]->c->ny,faces[i]->c->nz);
			glTexCoord2f((atan2(faces[i]->c->z +.25, faces[i]->c->x) + PI)/(2*PI), (faces[i]->c->y + 1.2)/2.4); 
			glVertex3f(faces[i]->c->x,faces[i]->c->y, faces[i]->c->z); 

			glNormal3f(faces[i]->d->nx,faces[i]->d->ny,faces[i]->d->nz);
			glTexCoord2f((atan2(faces[i]->d->z +.25, faces[i]->d->x) + PI)/(2*PI), (faces[i]->d->y + 1.2)/2.4); 
			glVertex3f(faces[i]->d->x,faces[i]->d->y, faces[i]->d->z); 
				
			glEnd();

		}

	}

	//swap buffers so no flickering
	glutSwapBuffers();
	glFlush ();
	
	glutPostRedisplay();

}
//keyboard function to interact with the program
void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
   		//togle wireframe
		case 'w':
			wireframe = !wireframe;
		break;
		//clear the mesh	
		case 'r':
			clearmesh();
			calcnormals();
		break;
		//subdivide mesh	
		case 'c':
			//code to subdivide
			verts = vertices;
			calcfaces();
			calcedges();
			calcverts();
			createmesh(); 
			calcnormals();
			//clear the previous mesh
			faces.clear();
			edges.clear();
			//set the current mesh to the subdived mesh	
			faces = tempfaces;
			edges = tempedges;
		break;
		//Escape Key - Exit Program	
		case 27:
			exit(0);
		break;
   }
}
//reshape function for when window is resized
void reshape (int w, int h)
{
	glViewport (0, 0,(GLsizei)w,(GLsizei)h); 
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluPerspective(90.0,(float)w/h,0.01,10.0);
	glMatrixMode (GL_MODELVIEW);
}

//main method
int main(int argc, char** argv)
{
	//Set Up Glut
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	//set Window size
	glutInitWindowSize (640, 480); 
	//Set window Position
	glutInitWindowPosition (50, 50);
	//Set window title
	glutCreateWindow ("Catmull Clark Subdivisions");
	//Set up the keyboard commands
	glutKeyboardFunc(keyboard);
	//set the display function
	glutDisplayFunc(display); 
	//Set the reshape function, for if the window is resized
	glutReshapeFunc(reshape);
 	//call the init method
	init();
	//start main loop
	glutMainLoop();
	return 0;
}