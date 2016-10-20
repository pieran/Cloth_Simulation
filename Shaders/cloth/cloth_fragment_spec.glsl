#version 150 core

uniform sampler2D diffuseTex;

uniform vec3  ambientColour;
uniform vec3  invLightDir;			//Directional Light
uniform vec3  cameraPos;
uniform float specularIntensity;

uniform float isTextured;
uniform float minVal;
uniform float maxVal;
uniform sampler1D colorGradient;


in Vertex	{
	vec3 worldPos;
	vec2 texCoord;
	float val;
	vec3 normal;
} IN;

out vec4 gl_FragColor;

void main(void)	{	
	

	   
//Colour Computations
	vec3 col;
	if (isTextured == 1.0f)
	{
		//Lighting Calculations
		vec3 normal 		= normalize(IN.normal);
		vec3 viewDir 		= normalize(cameraPos - IN.worldPos );
		vec3 halfDir 		= normalize(invLightDir + viewDir );
		float rFactor       = max(0.0, dot(halfDir , normal ));
		
		float dFactor       = max(0.0, dot(invLightDir , normal )) ;
		float sFactor       = pow(rFactor , specularIntensity );
	
		vec4 texColourA  = texture(diffuseTex, IN.texCoord * 5.0f);
		vec4 texColourB  = texture(diffuseTex, IN.texCoord);
		
		vec4 texColour  = mix(texColourA, texColourB, 0.4f);
		vec3 specColour = min(texColour.rgb + vec3(0.5f), vec3(1)); //Quick hack to approximate specular colour of an object, assuming the light colour is white
		col = ambientColour * texColour .rgb
				 + texColour.rgb * dFactor
				 + specColour * sFactor * 0.33f;
				 
		gl_FragColor.a 		= texColour.a;
	}
	else
	{
		float val 	= clamp((IN.val - minVal) / (maxVal - minVal), 0.0f, 1.0f);
		col 		= texture(colorGradient, val).rgb;					 
		gl_FragColor.a 		= 1.0f;
	}
		
//Output Final Fragment Colour
	gl_FragColor.rgb 	= col;
}