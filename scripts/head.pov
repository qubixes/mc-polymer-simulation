#include "shapes.inc"
#include "colors.inc"

#declare Saw_Node = sphere {
	0, 0.1
}

#declare Saw_Edge = cylinder {
	0, x, 0.05
}

global_settings {
	radiosity {
		pretrace_start 0.08
		pretrace_end 0.04
		count 35
		nearest_count 5
		error_bound 1.8
		recursion_limit 3
		low_error_factor 0.5
		gray_threshold 0.0
		minimum_reuse 0.015
		brightness 1
		
		adc_bailout 0.01/2
		}
	}

#declare Saw_Piece = merge {
	object { Saw_Edge }
	object { Saw_Node translate x }
}

// object { Saw_Node texture {pigment {Orange} } scale 2 }

// camera {
// 	location <31, 25, -59>
// 	focal_point 0
// 	look_at x-0.15*y
// 	angle 30
// 	aperture 0.2
// 	blur_samples 50
// }

// background { color Gray10 }
sky_sphere {
	pigment { color White }
}

// light_source { <29, 27, -9> color White area_light <10,0,0>, <0,0,10>, 3, 3
// 	adaptive 5
// 	jitter
// }

// light_source { <-43, 30, 22> color White area_light <5,0,0>, <0,0,5>, 5, 5
// 	adaptive 5
// // 	jitter
// }
