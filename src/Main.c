#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/Random.h"

#include "./Math3D.h"


mesh meshCube;

mat4x4 matProj;
vec3d vCamera = { 0.0f,0.0f,0.0f,1.0f };
vec3d vVelocity = { 0.0f,0.0f,0.0f,1.0f };
vec3d vLookDir;
float fYaw;	
float fPitch;	
float fTheta;

vec3d vLength = { 0.5f,1.8f,0.5f,1.0f };
Vec2 MouseBefore = { 0.0f,0.0f };

mat4x4 matView;

int Mode = 0;
int Menu = 0;

Vec2 FunctionOrigin = { 0.0f,0.0f };


float Function_2D(float x,float y){
	return (1.0f / (x * x + 1)) + (1.0f / (y * y + 1));
}

int compare(const void* e1,const void* e2) {
	triangle t1 = *(triangle*)e1;
	triangle t2 = *(triangle*)e2;
	float z1 = (t1.p[0].z+t1.p[1].z+t1.p[2].z)/3;
    float z2 = (t2.p[0].z+t2.p[1].z+t2.p[2].z)/3;
    return z1 == z2 ? 0 : (z1 < z2 ? 1 : -1);
}
void Triangles_Project(){
	Vector vecTrianglesToRaster = Vector_New(sizeof(triangle));

	for (int i = 0;i<meshCube.tris.size;i++){
		triangle tri = *(triangle*)Vector_Get(&meshCube.tris,i);
		
		vec3d vCameraRay = vec3d_Sub(tri.p[0], vCamera);

		if (vec3d_DotProduct(tri.n,vCameraRay) < 0.0f){
			tri.p[0] = Matrix_MultiplyVector(matView,tri.p[0]);
			tri.p[1] = Matrix_MultiplyVector(matView,tri.p[1]);
			tri.p[2] = Matrix_MultiplyVector(matView,tri.p[2]);

			int nClippedTriangles = 0;
			triangle clipped[2];
			nClippedTriangles = Triangle_ClipAgainstPlane(vec3d_new(0.0f,0.0f,0.1f),vec3d_new(0.0f,0.0f,1.0f),tri,&clipped[0],&clipped[1]);

			for (int n = 0; n < nClippedTriangles; n++){
				clipped[n].p[0] = Matrix_MultiplyVector(matProj, clipped[n].p[0]);
				clipped[n].p[1] = Matrix_MultiplyVector(matProj, clipped[n].p[1]);
				clipped[n].p[2] = Matrix_MultiplyVector(matProj, clipped[n].p[2]);

				clipped[n].p[0] = vec3d_Div(clipped[n].p[0], clipped[n].p[0].w);
				clipped[n].p[1] = vec3d_Div(clipped[n].p[1], clipped[n].p[1].w);
				clipped[n].p[2] = vec3d_Div(clipped[n].p[2], clipped[n].p[2].w);

				clipped[n].p[0].x *= -1.0f;
				clipped[n].p[1].x *= -1.0f;
				clipped[n].p[2].x *= -1.0f;
				clipped[n].p[0].y *= -1.0f;
				clipped[n].p[1].y *= -1.0f;
				clipped[n].p[2].y *= -1.0f;

				vec3d vOffsetView = vec3d_new( 1,1,0 );
				clipped[n].p[0] = vec3d_Add(clipped[n].p[0], vOffsetView);
				clipped[n].p[1] = vec3d_Add(clipped[n].p[1], vOffsetView);
				clipped[n].p[2] = vec3d_Add(clipped[n].p[2], vOffsetView);
				clipped[n].p[0].x *= 0.5f * (float)GetWidth();
				clipped[n].p[0].y *= 0.5f * (float)GetHeight();
				clipped[n].p[1].x *= 0.5f * (float)GetWidth();
				clipped[n].p[1].y *= 0.5f * (float)GetHeight();
				clipped[n].p[2].x *= 0.5f * (float)GetWidth();
				clipped[n].p[2].y *= 0.5f * (float)GetHeight();

				Vector_Push(&vecTrianglesToRaster,&clipped[n]);
			}			
		}
	}

	qsort(vecTrianglesToRaster.Memory,vecTrianglesToRaster.size,vecTrianglesToRaster.ELEMENT_SIZE,compare);

	for (int i = 0;i<vecTrianglesToRaster.size;i++)
	{
		triangle triToRaster = *(triangle*)Vector_Get(&vecTrianglesToRaster,i);

		triangle clipped[2];
		Vector listTriangles = Vector_New(sizeof(triangle));

		Vector_Push(&listTriangles,&triToRaster);
		int nNewTriangles = 1;

		for (int p = 0; p < 4; p++)
		{
			int nTrisToAdd = 0;
			while (nNewTriangles > 0)
			{
				triangle test = *(triangle*)Vector_Get(&listTriangles,0);
				Vector_Remove(&listTriangles,0);
				nNewTriangles--;

				switch (p)
				{
				case 0:	nTrisToAdd = Triangle_ClipAgainstPlane(vec3d_new( 0.0f, 0.0f, 0.0f ), 					vec3d_new( 0.0f, 1.0f, 0.0f ), 	test, &clipped[0], &clipped[1]); break;
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane(vec3d_new( 0.0f, (float)GetHeight() - 1, 0.0f ), vec3d_new( 0.0f, -1.0f, 0.0f ), test, &clipped[0], &clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane(vec3d_new( 0.0f, 0.0f, 0.0f ), 					vec3d_new( 1.0f, 0.0f, 0.0f ), 	test, &clipped[0], &clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane(vec3d_new( (float)GetWidth() - 1, 0.0f, 0.0f ), 	vec3d_new( -1.0f, 0.0f, 0.0f ), test, &clipped[0], &clipped[1]); break;
				}

				for (int w = 0; w < nTrisToAdd; w++)
					Vector_Push(&listTriangles,&clipped[w]);
			}
			nNewTriangles = listTriangles.size;
		}

		for (int j = 0;j<listTriangles.size;j++){
			triangle t = *(triangle*)Vector_Get(&listTriangles,j);
			//RenderTriangle(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c);
			//RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),WHITE,1.0f);

			if(Mode==0) RenderTriangle(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c);
			if(Mode==1) RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c,1.0f);
			if(Mode==2){
				RenderTriangle(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c);
				RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),WHITE,1.0f);
			}
		}

		Vector_Free(&listTriangles);
	}
	Vector_Free(&vecTrianglesToRaster);
}
void Menu_Set(int m){
	if(Menu==0 && m==1){
		AlxWindow_Mouse_SetInvisible(&window);
		SetMouse((Vec2){ GetWidth() / 2,GetHeight() / 2 });
	}
	if(Menu==1 && m==0){
		AlxWindow_Mouse_SetVisible(&window);
	}
	
	MouseBefore = GetMouse();
	Menu = m;
}


void Setup(AlxWindow* w){
	Menu_Set(1);

	RGA_Set(Time_Nano());
	RGA_Get(6969);

	meshCube = (mesh){ Vector_New(sizeof(triangle)) };
	matProj = Matrix_MakeProjection(90.0f, (float)GetHeight() / (float)GetWidth(), 0.1f, 1000.0f);
}

void Update(AlxWindow* w){
	if(Menu==1){
		if(GetMouse().x!=MouseBefore.x || GetMouse().y!=MouseBefore.y){
			Vec2 d = Vec2_Sub(GetMouse(),MouseBefore);
			Vec2 a = Vec2_Mulf(Vec2_Div(d,(Vec2){ window.Width,window.Height }),2 * F64_PI);
	
			fYaw += a.x;
			fPitch += a.y;
	
			SetMouse((Vec2){ GetWidth() / 2,GetHeight() / 2 });
			MouseBefore = GetMouse();
		}
	}
	
	if(Stroke(ALX_KEY_ESC).PRESSED)
		Menu_Set(!Menu);

	if(Stroke(ALX_KEY_UP).DOWN)
		fPitch -= 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_DOWN).DOWN)
		fPitch += 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_LEFT).DOWN)
		fYaw -= 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_RIGHT).DOWN)
		fYaw += 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_Z).PRESSED)
		Mode = Mode < 3 ? Mode+1 : 0;

	
	if(Stroke(ALX_KEY_R).DOWN)
		vVelocity.y = 5.0f;
	else if(Stroke(ALX_KEY_F).DOWN)
		vVelocity.y = -5.0f;
	else
		vVelocity.y = 0.0f;

	mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
	vec3d vForward = Matrix_MultiplyVector(matCameraRot,vec3d_new(0.0f,0.0f,1.0f));
	vec3d vLeft = vec3d_Perp(vForward);
	
	if(Stroke(ALX_KEY_W).DOWN){
		vCamera = vec3d_Add(vCamera,vec3d_Mul(vForward,5.0f * w->ElapsedTime));
	}
	if(Stroke(ALX_KEY_S).DOWN){
		vCamera = vec3d_Add(vCamera,vec3d_Mul(vForward,-5.0f * w->ElapsedTime));
	}
	if(Stroke(ALX_KEY_A).DOWN){
		vCamera = vec3d_Add(vCamera,vec3d_Mul(vLeft,-5.0f * w->ElapsedTime));
	}
	if (Stroke(ALX_KEY_D).DOWN){
		vCamera = vec3d_Add(vCamera,vec3d_Mul(vLeft,5.0f * w->ElapsedTime));
	}

	vVelocity = vec3d_Add(vVelocity,vec3d_Mul((vec3d){ 0.0f,0.0f,0.0f,1.0f },w->ElapsedTime));
	vCamera = vec3d_Add(vCamera,vec3d_Mul(vVelocity,w->ElapsedTime));

	float Border = F32_PI * 0.5f - 0.001f;
	if(fPitch<-Border) fPitch = -Border;
	if(fPitch>Border) fPitch = Border;

	vec3d vUp = vec3d_new( 0.0f,1.0f,0.0f );
	vec3d vTarget = vec3d_new( 0.0f,0.0f,1.0f );
	mat4x4 matCameraRotX = Matrix_MakeRotationX(fPitch);
	vLookDir = Matrix_MultiplyVector(matCameraRotX,vTarget);
	vLookDir = Matrix_MultiplyVector(matCameraRot,vLookDir);
	
	vTarget = vec3d_Add(vCamera, vLookDir);
	mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);
	matView = Matrix_QuickInverse(matCamera);


	Vector_Clear(&meshCube.tris);
	for(float i = -7;i<7;i+=0.2f){
		for(float j = -7;j<7;j+=0.2f){
			const float x1 = FunctionOrigin.x + j;
			const float x2 = FunctionOrigin.x + j+1;
			const float y1 = FunctionOrigin.y + i;
			const float y2 = FunctionOrigin.y + i+1;

			const vec3d p1 = { x1,Function_2D(x1,y1),y1,1.0f };
			const vec3d p2 = { x2,Function_2D(x2,y1),y1,1.0f };
			const vec3d p3 = { x1,Function_2D(x1,y2),y2,1.0f };
			const vec3d p4 = { x2,Function_2D(x2,y2),y2,1.0f };

			triangle t1 = { p1,p2,p4,{},RED };
			triangle t2 = { p1,p4,p3,{},RED };
			triangle t3 = { p1,p4,p2,{},RED };
			triangle t4 = { p1,p3,p4,{},RED };

			Triangle_CalcNorm(&t1);
			Triangle_CalcNorm(&t2);
			Triangle_CalcNorm(&t3);
			Triangle_CalcNorm(&t4);

			Triangle_ShadeNorm(&t1,(vec3d){ 0.5f,0.4f,0.6f,1.0f });
			Triangle_ShadeNorm(&t2,(vec3d){ 0.5f,0.4f,0.6f,1.0f });
			Triangle_ShadeNorm(&t3,(vec3d){ 0.5f,0.4f,0.6f,1.0f });
			Triangle_ShadeNorm(&t4,(vec3d){ 0.5f,0.4f,0.6f,1.0f });

			Vector_Push(&meshCube.tris,&t1);
			Vector_Push(&meshCube.tris,&t2);
			Vector_Push(&meshCube.tris,&t3);
			Vector_Push(&meshCube.tris,&t4);
		}
	}


	Clear(LIGHT_BLUE);
	Triangles_Project();

	String str = String_Format("X: %f, Y: %f, Z: %f, Size: %d",vCamera.x,vCamera.y,vCamera.z,meshCube.tris.size);
	RenderCStrSize(str.Memory,str.size,0,0,RED);
	String_Free(&str);
}

void Delete(AlxWindow* w){
	Vector_Free(&meshCube.tris);
	AlxWindow_Mouse_SetVisible(&window);
}

int main(){
    if(Create("Geogebra 3D 1.0",2500,1200,1,1,Setup,Update,Delete))
        Start();
    return 0;
}