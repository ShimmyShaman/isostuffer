package isostuffer

import "core:fmt"
import mx "core:math"
import la "core:math/linalg"
import q "core:container/queue"

vec3 :: la.Vector3f32

ProcSDF :: #type proc(pos: vec3) -> f32

Point3 :: struct {
  x, y, z: int,
}

TetrahedralMesh :: struct {

}

generate_tetrahedral_mesh :: proc(sdf: ProcSDF) -> (mesh: ^TetrahedralMesh, okay: bool = true) {
  // Create Boundary Octants
//     case 0: // create boundary octants
//     if ( !dvproc_ ) return;

//     delete pstuffer_;
//     dvproc_->update_mesh_params(params_.resolution(), 
//             params_.num_levels(), params_.margin());
//     canvas->scale_ = dvproc_->grid_size() * 0.15;
//     pstuffer_ = new TIsoStuffer(dvproc_->resolution(),
//             params_.num_levels(), dvproc_->grid_size(),
//             dvproc_->min_point(), params_.alpha_long(),
//             params_.alpha_short(), *pdistfunc_);

//     pstuffer_->create_boundary_octants();
//     ++ curStep_;
//     canvas->updateGL();
//     return;
  _create_boundary_octants(sdf)

  return nil, false

// case 1: // create octree
//     pstuffer_->create_octree();
//     std::cerr << "CHECK CHILDREN MASK ... I" << std::endl;
//     pstuffer_->oct_tree().check_children_mask();

//     pstuffer_->fill_highest_level();
//     //std::cerr << "CHECK CHILDREN MASK ... II" << std::endl;
//     //pstuffer_->oct_tree().check_children_mask();
//     ++ curStep_;
//     canvas->updateGL();
//     return;
// case 2: // weak balance
//     pstuffer_->weak_balance();
//     std::cerr << "CHECK CHILDREN MASK ... III" << std::endl;
//     pstuffer_->oct_tree().check_children_mask();
//     pstuffer_->update_bitmasks();
//     ++ curStep_;
//     canvas->updateGL();
//     return;
// case 3:
//     pstuffer_->create_background_tets();

//     //// create background mesh
//     delete pmesh_;
//     pmesh_ = new TetMesh<double>;
//     pstuffer_->background_tet_mesh(pmesh_);
//     pmesh_->update_surface();

//     ++ curStep_;
//     canvas->updateGL();
//     return;
// case 4:
//     pstuffer_->compute_cutting_pts();
//     pstuffer_->wrap_violated_pts();
//     ++ curStep_;
//     canvas->updateGL();
//     return;
// case 5:
//     printf("INFO: extract final tets...\n");
//     pstuffer_->extract_final_tets();

//     delete pmesh_;
//     pmesh_ = new TetMesh<double>;
//     pstuffer_->extract_final_mesh(pmesh_);
//     pmesh_->update_surface();

//     ++ curStep_;
//     canvas->updateGL();
//     return;
// case 6:
//     return;
// default:
//     fprintf(stderr, "Unknown step ID %d\n", curStep_);
//     SHOULD_NEVER_HAPPEN(1);
// }
}

LEAF_OCTANT_WIDTH: f32 = 1.0
LEAF_OCTANT_HALF_WIDTH: f32 = LEAF_OCTANT_WIDTH * 0.5

_find_initial_surface_intersecting_octant :: proc(sdf: ProcSDF) -> Point3 {
  non_negative, non_positive := false, false

  pt: Point3 = Point3 {0, 0, 0}
  initial: for pt.x = -32; pt.x < 32; pt.x += 1 {
    for pt.y = -32; pt.y < 32; pt.y += 1 {
      for pt.z = -32; pt.z < 32; pt.z += 1 {
        // Centre Point
        spc := mx.sign_f32(sdf(vec3{ auto_cast pt.x, auto_cast pt.y, auto_cast pt.z }))
        // fmt.println("spc:", spc)
        if spc >= 0 do non_negative = true
        if spc <= 0 do non_positive = true
    
        if non_negative && non_positive do return pt
    
        // Corner Points
        for x := 0; x < 2; x += 1 {
          for y := 0; y < 2; y += 1 {
            for z := 0; z < 2; z += 1 {
              sv := mx.sign_f32(sdf(vec3{ auto_cast pt.x - LEAF_OCTANT_HALF_WIDTH + auto_cast x * LEAF_OCTANT_WIDTH,
                                          auto_cast pt.y - LEAF_OCTANT_HALF_WIDTH + auto_cast y * LEAF_OCTANT_WIDTH,
                                          auto_cast pt.z - LEAF_OCTANT_HALF_WIDTH + auto_cast z * LEAF_OCTANT_WIDTH }))
              // fmt.println("sv:", sv)
                                          
              if sv >= 0 do non_negative = true
              if sv <= 0 do non_positive = true
              if non_negative && non_positive do return pt
            }
          }
        }
      }
    }
  }

  return Point3 {10000, 10000, 10000}
}

@(private) _bfs_boundary_octants :: proc(sdf: ProcSDF, initial_boundary_intersection: Point3) {
  // Create the collections
  surf_con: q.Queue(Point3)
  q.init(&surf_con)
  defer q.destroy(&surf_con)
  
  surf_checked := hashset_create()
  defer hashset_destroy(&surf_checked)
  djb2 :: proc(item: ^Point3) -> (result: int) {
    result = 5381
    result = ((result << 5) + result) + auto_cast item.x
    result = ((result << 5) + result) + auto_cast item.y
    result = ((result << 5) + result) + auto_cast item.z
    
    return result
  }
  
  non_negative, non_positive := false, false
  
  // Begin the queue
  q.push_back(&surf_con, initial_boundary_intersection)
  octree_create_leaf(&octree, initial_boundary_intersection)
  hashset_add(&surf_checked, djb2(&initial_boundary_intersection))

  // Iterate
  ok: bool
  pt: Point3
  for {
    pt, ok = q.pop_front_safe(&surf_con)
    if !ok do break

    // Centre Point
    spc := mx.sign_f32(sdf(vec3{ auto_cast pt.x, auto_cast pt.y, auto_cast pt.z }))
    if spc >= 0 do non_negative = true
    if spc <= 0 do non_positive = true

    /*
    * If a leaf octant o has a corner vertex v whose sign
    * is opposite the sign of o’s center point, or if either sign is
    * zero, then we must create the three leaf octants incident
    * on v that share a square face with o
    */
  
    // Corner Points
    scc: [8]f32
    for x := 0; x < 2; x += 1 {
      for y := 0; y < 2; y += 1 {
        for z := 0; z < 2; z += 1 {
          sv := mx.sign_f32(sdf(vec3{  auto_cast pt.x - LEAF_OCTANT_HALF_WIDTH + auto_cast x * LEAF_OCTANT_WIDTH,
                                      auto_cast pt.y - LEAF_OCTANT_HALF_WIDTH + auto_cast y * LEAF_OCTANT_WIDTH,
                                      auto_cast pt.z - LEAF_OCTANT_HALF_WIDTH + auto_cast z * LEAF_OCTANT_WIDTH }))
          scc[x + y * 2 + z * 4] = sv

          if spc * sv <= 0 {
            //         // NOTE: here nowpos.x+D_COVTX[i][0], nowpos.y+D_.., nowpos.z+D_.. should NOT be out of range
            //         // create three octants 
            //         octTree_.create_leaf(nowpos.x+D_COVTX[i][0], nowpos.y, nowpos.z);
            //         octTree_.create_leaf(nowpos.x, nowpos.y+D_COVTX[i][1], nowpos.z);
            //         octTree_.create_leaf(nowpos.x, nowpos.y, nowpos.z+D_COVTX[i][2]);
            // Create
          }
        }
      }
    }
    //     if ( isovalue * ivs[i] <= 0. ) 
    //     {   
    //         // NOTE: here nowpos.x+D_COVTX[i][0], nowpos.y+D_.., nowpos.z+D_.. should NOT be out of range
    //         // create three octants 
    //         octTree_.create_leaf(nowpos.x+D_COVTX[i][0], nowpos.y, nowpos.z);
    //         octTree_.create_leaf(nowpos.x, nowpos.y+D_COVTX[i][1], nowpos.z);
    //         octTree_.create_leaf(nowpos.x, nowpos.y, nowpos.z+D_COVTX[i][2]);
    //     }
    // }

    // for(int i = 0;i < 12;++ i)
    //     cross[i] = ivs[IDX_EDGES[i][0]] * ivs[IDX_EDGES[i][1]] < 0;

    // for(int i = 0;i < 3;++ i)
    // for(int j = 0, dir = -1;j < 2;++ j, dir *= -1)
    // {
    //     Tuple3i newp = nowpos;
    //     newp[i] += dir;

    //     if ( newp[i] >= 0 && newp[i] < highestRes_[i] && !visited[newp.z][newp.y][newp.x] )
    //     {
    //         const int fid = (i<<1) + j;
    //         int k;
    //         for(k = 0;k < 4;++ k) if ( cross[IDX_EDGE_FACES[fid][k]] ) break;

    //         if ( k < 4 )
    //         {
    //             visited[newp.z][newp.y][newp.x] = true;
    //             octTree_.create_leaf(newp.x, newp.y, newp.z);
    //             que.push(newp);
    //         }
    //     }
    // }

    if !non_negative || !non_positive do continue

    // Add the unit octants
    for i := 0; i < 8; i += 1 {
      if spc * scc[i] <= 0.0 {
        // NOTE: here nowpos.x+D_COVTX[i][0], nowpos.y+D_.., nowpos.z+D_.. should NOT be out of range
        // create three octants 
        octTree_.create_leaf(nowpos.x+D_COVTX[i][0], nowpos.y, nowpos.z);
        octTree_.create_leaf(nowpos.x, nowpos.y+D_COVTX[i][1], nowpos.z);
        octTree_.create_leaf(nowpos.x, nowpos.y, nowpos.z+D_COVTX[i][2]);
      }
    }
//         for(int i = 0;i < 12;++ i)
//             cross[i] = ivs[IDX_EDGES[i][0]] * ivs[IDX_EDGES[i][1]] < 0;

//         for(int i = 0;i < 3;++ i)
//         for(int j = 0, dir = -1;j < 2;++ j, dir *= -1)
//         {
//             Tuple3i newp = nowpos;
//             newp[i] += dir;

//             if ( newp[i] >= 0 && newp[i] < highestRes_[i] && !visited[newp.z][newp.y][newp.x] )
//             {
//                 const int fid = (i<<1) + j;
//                 int k;
//                 for(k = 0;k < 4;++ k) if ( cross[IDX_EDGE_FACES[fid][k]] ) break;

//                 if ( k < 4 )
//                 {
//                     visited[newp.z][newp.y][newp.x] = true;
//                     octTree_.create_leaf(newp.x, newp.y, newp.z);
//                     que.push(newp);
//                 }
//             }
//         }
  }






/*
 * if a leaf octant o has a corner vertex v whose sign
 * is opposite the sign of o’s center point, or if either sign is
 * zero, then we must create the three leaf octants incident
 * on v that share a square face with o
 */
// template <class DistFunc>
// void IsoStuffer<DistFunc>::bfs_boundary_octants(int ix, int iy, int iz)
// {
//     using namespace std;

//     queue<Tuple3i> que;
//     TLatticeb visited(boost::extents[highestRes_.z][highestRes_.y][highestRes_.x]);
//     zero_multi_array(visited);
    
//     que.push(Tuple3i(ix, iy, iz));
//     visited[iz][iy][ix] = true;
//     octTree_.create_leaf(ix, iy, iz);

//     double ivs[8];
//     bool   cross[12];   // indicate if there is a cutting point on each of 12 edges 

//     while ( !que.empty() )
//     {
//         Tuple3i nowpos = que.front();
//         que.pop();

//         double isovalue = get_center_iso_value(nowpos);

//         for(int i = 0;i < 8;++ i)
//         {
//             ivs[i]= get_iso_value(nowpos.x+D_CORNERS[i][0], 
//                                   nowpos.y+D_CORNERS[i][1], 
//                                   nowpos.z+D_CORNERS[i][2]);

//             if ( isovalue * ivs[i] <= 0. ) 
//             {   
//                 // NOTE: here nowpos.x+D_COVTX[i][0], nowpos.y+D_.., nowpos.z+D_.. should NOT be out of range
//                 // create three octants 
//                 octTree_.create_leaf(nowpos.x+D_COVTX[i][0], nowpos.y, nowpos.z);
//                 octTree_.create_leaf(nowpos.x, nowpos.y+D_COVTX[i][1], nowpos.z);
//                 octTree_.create_leaf(nowpos.x, nowpos.y, nowpos.z+D_COVTX[i][2]);
//             }
//         }

//         for(int i = 0;i < 12;++ i)
//             cross[i] = ivs[IDX_EDGES[i][0]] * ivs[IDX_EDGES[i][1]] < 0;

//         for(int i = 0;i < 3;++ i)
//         for(int j = 0, dir = -1;j < 2;++ j, dir *= -1)
//         {
//             Tuple3i newp = nowpos;
//             newp[i] += dir;

//             if ( newp[i] >= 0 && newp[i] < highestRes_[i] && !visited[newp.z][newp.y][newp.x] )
//             {
//                 const int fid = (i<<1) + j;
//                 int k;
//                 for(k = 0;k < 4;++ k) if ( cross[IDX_EDGE_FACES[fid][k]] ) break;

//                 if ( k < 4 )
//                 {
//                     visited[newp.z][newp.y][newp.x] = true;
//                     octTree_.create_leaf(newp.x, newp.y, newp.z);
//                     que.push(newp);
//                 }
//             }
//         }
//     } // end while 
// }
}

@(private) _create_boundary_octants :: proc(sdf: ProcSDF) -> (okay: bool = true) {
  fmt.println("create boundary octants")

  // surf_checked := hashset_create()
  // defer hashset_destroy(&surf_checked)
  // djb2 :: proc(item: ^Point3) -> (result: int) {
  //   result = 5381
  //   result = ((result << 5) + result) + auto_cast item.x
  //   result = ((result << 5) + result) + auto_cast item.y
  //   result = ((result << 5) + result) + auto_cast item.z
    
  //   return result
  // }
  // surf_con: q.Queue(Point3)
  // q.init(&surf_con)
  // defer q.destroy(&surf_con)

  // -- Find a surface-intersecting octant by sweeping the grid, expanding from the origin
  pt: Point3 = _find_initial_surface_intersecting_octant(sdf)
  if pt.x == 10000 {
    fmt.println("ERROR: Could not find a surface-intersecting octant")
    okay = false
    return
  }

  _bfs_boundary_octants(sdf, pt)

  fmt.println("initial surface intersecting octant:", pt)
  return
  
  // sio: for {
  //   added, okay := _check_octant(&octree, sdf, pt)
  //   if !okay {
  //     fmt.println("ERROR: Checking initial surface-intersecting octant")
  //     return
  //   }
  //   if added {
  //     // Add the first surface point
  //     break
  //   }

  //   // -- Move surf_con from the origin
  //   pt = Point3 {pt.x + 1, pt.y + 1, pt.z + 1}
  //   if pt.x >= 32 {
  //     fmt.println("ERROR: Could not find a surface-intersecting octant")
  //     okay = false
  //     return
  //   }
  // }
  
  // // Centre Point
  // spc := mx.sign_f32(sdf(vec3{ auto_cast pt.x, auto_cast pt.y, auto_cast pt.z }))
  // if spc >= 0 do non_negative = true
  // if spc <= 0 do non_positive = true

  // // Corner Points
  // scc: [8]f32
  // for x := 0; x < 2; x += 1 {
  //   for y := 0; y < 2; y += 1 {
  //     for z := 0; z < 2; z += 1 {
  //       sv := mx.sign_f32(sdf(vec3{  auto_cast pt.x - LEAF_OCTANT_HALF_WIDTH + auto_cast x * LEAF_OCTANT_WIDTH,
  //                                   auto_cast pt.y - LEAF_OCTANT_HALF_WIDTH + auto_cast y * LEAF_OCTANT_WIDTH,
  //                                   auto_cast pt.z - LEAF_OCTANT_HALF_WIDTH + auto_cast z * LEAF_OCTANT_WIDTH }))

  //       if sv >= 0 do non_negative = true
  //       if sv <= 0 do non_positive = true
        
  //       scc[x + y * 2 + z * 4] = sv
  //     }
  //   }
  // }

  // if !non_negative || !non_positive do return

  // added = true
  // if add_to_octree(octree, pt, 1) {
  //   okay = false
  //   return
  // }


  // hashset_add(&surf_checked, djb2(&pt))
  // q.push_back(&surf_con, pt)
}