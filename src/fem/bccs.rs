use std::mem;
use crate::fem::mesh::Mesh;
use crate::fem::msg::Messenger;
use crate::fem::sparse::SparseMatrix;
use crate::fem::solver::GenericSolver;
use super::error::FemError;


pub struct BccSparseMatrix {
    num_vertex: usize,
    block_size: usize,
    ptr:        Vec<usize>,
    ind:        Vec<usize>,
    val:        Vec<f64>,
}

struct BccFactor {
    num_vertex: usize,
    block_size: usize,
    l_space:    usize,
    l_ind:      Vec<usize>,
    perm:       Vec<usize>,
    inv_p:      Vec<usize>,
    x_ind:      Vec<usize>,
    x_val:      Vec<usize>,
    d_val:      Vec<f64>,
    l_val:      Vec<f64>,
}
impl BccFactor {
    pub fn new(num_vertex: usize, block_size: usize) -> Self {
        Self {
            num_vertex,
            block_size,
            l_space: 0,
            l_ind: vec![],
            perm: vec![],
            inv_p: vec![],
            x_ind: vec![],
            x_val: vec![],
            d_val: vec![],
            l_val: vec![]
        }
    }
}

impl BccSparseMatrix {
    pub fn new(mesh: &Mesh) -> Result<Self, FemError> {
        let mut ptr: Vec<usize> = vec!();
        let mut ind: Vec<usize> = vec!();
        let mut val: Vec<f64> = vec!();
        Self::set_matrix(mesh, &mut ptr, &mut ind, &mut val);
        Ok(Self { num_vertex: mesh.num_vertex, block_size: mesh.freedom, ptr, ind, val })
    }
    fn set_matrix(mesh: &Mesh, ptr: &mut Vec<usize>, ind: &mut Vec<usize>, val: &mut Vec<f64>) {
        // Allocate & set all element adjacency lengths to null
        let mut n_ptr: Vec<usize> = vec![0; mesh.num_vertex + 1];

        // Calculate sizes for a list of all references to adjacency elements
        for i in 0..mesh.fe.shape()[0] * mesh.fe.shape()[1] {
            n_ptr[mesh.fe[[i / mesh.fe.shape()[1], i % mesh.fe.shape()[1]]]] += 1;
        }

        // Calculate references from lengths
        let mut node = 0;
        for i in 0..mesh.num_vertex {
            let j = n_ptr[i];
            n_ptr[i] = node;
            node += j;
        }
        n_ptr[mesh.num_vertex] = node;
        let mut n_ind = vec![0; node];

        // Save into n_ind all references to adjacency elements
        let mut k = 0;
        for i in 0..mesh.fe.shape()[0] {
            for _ in 0..mesh.fe.shape()[1] {
                n_ind[n_ptr[mesh.fe[[k / mesh.fe.shape()[1], k % mesh.fe.shape()[1]]]]] = i;
                n_ptr[mesh.fe[[k / mesh.fe.shape()[1], k % mesh.fe.shape()[1]]]] += 1;
                k += 1;
            }
        }

        // Restore n_ptr
        for i in (1..mesh.num_vertex + 1).rev() {
            n_ptr[i] = n_ptr[i - 1];
        }
        n_ptr[0] = 0;

        let mut mark:Vec<Option<usize>> = vec![None; mesh.num_vertex];
        ptr.resize(mesh.num_vertex + 1, 0);

        // Calculate memory for adjacency structure
        let mut n_edges = 0;
        for i in 0..mesh.num_vertex {
            ptr[i] = n_edges;
            node = n_ptr[i + 1];
            for j in n_ptr[i]..node {
                let mut m = mesh.fe.shape()[1] * n_ind[j];
                for _ in 0..mesh.fe.shape()[1] {
                    let n = mesh.fe[[m / mesh.fe.shape()[1], m % mesh.fe.shape()[1]]];
                    m += 1;
                    if mark[n] != Some(i) {
                        mark[n] = Some(i);
                        n_edges += 1;
                    }
                }
            }
        }
        ptr[mesh.num_vertex] = n_edges;
        ind.resize(n_edges, 0);

        // Save adjacency structure
        n_edges = 0;
        for i in 0..mesh.num_vertex {
            node = n_ptr[i + 1];
            for j in n_ptr[i]..node {
                let mut m = mesh.fe.shape()[1] * n_ind[j];
                for _ in 0..mesh.fe.shape()[1] {
                    let n = mesh.fe[[m / mesh.fe.shape()[1], m % mesh.fe.shape()[1]]];
                    m += 1;
                    if mark[n] != Some(i) {
                        mark[n] = Some(i);
                        ind[n_edges] = n;
                        n_edges += 1;
                    }
                }
            }
            (& mut ind[ptr[i]..n_edges]).sort();
        }
        val.resize(ptr[mesh.num_vertex] * mesh.freedom * mesh.freedom, 0.0);
    }
}
impl SparseMatrix for BccSparseMatrix {
    fn add_value(&mut self, i: usize, j: usize, value: f64) -> Result<(), FemError> {
        let offset = self.block_size * (i % self.block_size) + (j % self.block_size);
        let start = self.ptr[i / self.block_size];
        let stop = self.ptr[i / self.block_size+1];
        for k in start..stop {
            if self.ind[k] == j / self.block_size {
                self.val[self.block_size * self.block_size * k + offset] += value;
                break
            }
        }
        Ok(())
    }

    fn solve(&mut self, rhs: &Vec<f64>, tolerance: f64) -> Result<Vec<f64>, FemError> {
        let mut factor = BccFactor::new(self.num_vertex, self.block_size);
        sp_order(&mut factor, self)?;
        sp_factor(&mut factor, self, tolerance)?;
        let mut result = rhs.clone();
        sp_solve(&mut factor, & mut result);
        Ok(result)
    }
    fn clear(&mut self) {
        self.ptr.clear();
        self.ind.clear();
        self.val.clear();
    }
    fn reset(&mut self) {
        // self.ptr.clear();
        // self.ind.clear();
        self.val.fill(0.0);
    }
}

pub fn new_bcc_solver(mesh: &Mesh) -> Result<GenericSolver<BccSparseMatrix>, FemError> {
    Ok(GenericSolver::new(BccSparseMatrix::new(mesh)?, mesh.num_vertex * mesh.freedom))
}

fn sp_order(factor: &mut BccFactor, matrix: &BccSparseMatrix) -> Result<(), FemError> {
    // Allocate standard graph structure without diagonal indices
    let mut x_adj = vec![0; matrix.num_vertex + 1];
    let mut adj = vec![0; matrix.ptr[matrix.num_vertex] - matrix.num_vertex];

    let mut msg = Messenger::new("Preparing system of equation", 0, factor.num_vertex as i64, 5);
    let mut l = 0;
    for i in 0..matrix.num_vertex {
        // if *canceled {
        //     return fmt.Errorf("the process was stopped by user")
        // }
        msg.add_progress();

        x_adj[i] = l;
        let start = matrix.ptr[i];
        let stop = matrix.ptr[i + 1];
        for j in start..stop {
            let k = matrix.ind[j];
            if k != i {
                adj[l] = k;
                l += 1;
            }
        }
    }
    x_adj[matrix.num_vertex] = l;

    // Now compute permutation
    factor.perm.resize(matrix.num_vertex, 0);
    factor.inv_p.resize(matrix.num_vertex, 0);
    gen_qmd(matrix.num_vertex, &mut x_adj, &mut adj, &mut factor.perm, &mut factor.inv_p);

    // Perform symbolic factorization
    sym_factorize(factor, &matrix.ptr, &matrix.ind);
    Ok(())
}

fn gen_qmd(num_vertex: usize, x_adj: &mut Vec<usize>, adj: &mut Vec<usize>, perm: &mut Vec<usize>, inv: &mut Vec<usize>) {
    // Initialisation
    let mut num = 0_usize;
    let mut search = 0;
    let mut min_deg = num_vertex as i32;
    let mut q_head = vec![0; 7 * num_vertex];
    let mut q_size  = vec![0; 7 * num_vertex];
    let mut q_link  = vec![0; 7 * num_vertex];
    let mut marker  = vec![0; 7 * num_vertex];
    let mut degree :Vec<i32> = vec![0; 7 * num_vertex];
    let mut reach  = vec![0; 7 * num_vertex];
    let mut neighbor  = vec![0; 7 * num_vertex];

    for i in 0..num_vertex {
        perm[i] = i;
        inv[i] = i;
        q_head[i] = i;
        q_size[i] = 1;
        q_link[i] = usize::MAX;
        marker[i] = 0;
        let n_deg = (x_adj[i + 1] - x_adj[i]) as i32;
        degree[i] = n_deg;
        if n_deg < min_deg {
            min_deg = n_deg
        }
    }
    let mut thresh = min_deg;
    min_deg = num_vertex as i32;

    while num < num_vertex {
        // Threshold searching
        loop {
            if num >= search {
                search = num;
            }
            let mut flag = 1;
            for i in search..num_vertex {
                let node = perm[i];
                if marker[node] < 0 {
                    continue
                }
                let n_deg = degree[node];
                if n_deg <= thresh {
                    search = i;
                    flag = 0;
                    break;
                }
                if n_deg < min_deg {
                    min_deg = n_deg
                }
            }
            if flag == 0 {
                break
            }

            thresh = min_deg;
            min_deg = num_vertex as i32;
            search = 0
        }
        // Root is pretender!
        let root = perm[search];
        marker[root] = 1;
        let mut deg = degree[root];

        // Calc reachable power (set_power)
        let mut set_power = 0;
        let mut reach_size = 0;
        let mut neighbor_size = 0;
        let mut j_start = x_adj[root];
        let mut j_stop = x_adj[root + 1];
        if j_start < j_stop {
            for j in j_start..j_stop {
                let mut node = adj[j];
                if marker[node] != 0 {
                    continue
                }
                if degree[node] >= 0 {
                    reach[reach_size] = node;
                    reach_size += 1;
                    marker[node] = 1;
                    set_power += q_size[node];
                    continue
                }
                neighbor[neighbor_size] = node;
                neighbor_size += 1;
                marker[node] = -1;
                let mut k = x_adj[node];
                let mut k_stop = x_adj[node + 1];
                while k < k_stop {
                    node = adj[k];
                    k += 1;
                    if node >= num_vertex {
                        node -= num_vertex;
                        if node < num_vertex {
                            k = x_adj[node];
                            k_stop = x_adj[node + 1];
                            continue
                        }
                    break
                    }
                    if marker[node] != 0 {
                        continue
                    }
                    reach[reach_size] = node;
                    reach_size += 1;
                    set_power += q_size[node];
                    marker[node] = 1;
                }
            }
        }

        // Numerate merged roots
        let mut nx_node = root;
        while nx_node < num_vertex {
            let np = inv[nx_node];
            let ip = perm[num];
            perm[np] = ip;
            inv[ip] = np;
            perm[num] = nx_node;
            inv[nx_node] = num;
            degree[nx_node] = -1;
            deg -= 1;
            num += 1;
            nx_node = q_link[nx_node];
        }

        if num >= num_vertex {
        break
        }
        if reach_size == 0 {
            // Isolate root
            marker[root] = 0;
            continue
        }
        let mut counter = 0;

        // Find merged roots
        let mut sn_size = 0;
        for i in 0..reach_size {
            let node = reach[i];
            j_start = x_adj[node];
            j_stop = x_adj[node + 1];
            for j in j_start..j_stop {
                let index = adj[j];
                if degree[index] >= 0 {
                    continue
                }
                if marker[index] != 0 {
                    continue
                }
                neighbor[neighbor_size + sn_size] = index;
                sn_size += 1;
                marker[index] = -1;
            }
        }
        if sn_size != 0 {
            for i in 0..sn_size {
                marker[neighbor[neighbor_size + i]] = 0;
            }
            for i in 0..sn_size {
                let tmp_root= neighbor[neighbor_size + i];
                marker[tmp_root] = -1;
                deg = 0;
                let mut ovl_size = 0;
                let mut search_size = 0;
                let mut node = tmp_root;
                let mut j = x_adj[node];
                j_stop = x_adj[node + 1];
                while j < j_stop {
                    node = adj[j];
                    j += 1;
                    if node >= num_vertex {
                        node -= num_vertex;
                        if node < num_vertex {
                            j = x_adj[node];
                            j_stop = x_adj[node + 1];
                            continue
                        }
                        break
                    }
                    if marker[node] < 0 {
                        continue
                    }
                    if marker[node] > 0 {
                        if marker[node] > 1 {
                            continue
                        }
                        marker[node] = 2; // Marker for overlapped
                        neighbor[neighbor_size + sn_size + ovl_size] = node;
                        ovl_size += 1;
                        continue
                    }
                    marker[node] = -2;
                    reach[reach_size + search_size] = node;
                    search_size += 1;
                    deg += q_size[node];
                }

                let mut head = num_vertex;
                let mut merge_size = 0;
                for j in 0..ovl_size {
                    node = neighbor[neighbor_size + sn_size + j];
                    let k_start = x_adj[node];
                    let k_stop = x_adj[node + 1];
                    let mut flag = 0;
                    for k in k_start..k_stop {
                        if marker[adj[k]] != 0 {
                            continue
                        }
                        flag = 1;
                        break
                    }
                    if flag != 0 {
                        marker[node] = 1;
                    } else {
                        merge_size += q_size[node];
                        counter += 1;
                        marker[node] = -1;
                        let l_node = q_head[node];
                        q_link[l_node] = head;
                        if head < num_vertex {
                            q_head[node] = q_head[head]
                        }
                        head = node;
                    }
                }
                if head < num_vertex {
                    degree[head] = set_power + deg - 1;
                    marker[head] = 2;
                    q_size[head] = merge_size;
                }
                marker[tmp_root] = 0;
                if search_size != 0 {
                    for j in 0..search_size {
                        marker[reach[reach_size+j]] = 0;
                    }
                }
            }
        }

        // Update root degree
        for i in 0..reach_size {
            let tmp_root = reach[i];
            if marker[tmp_root] > 1 || marker[tmp_root] < 0 {
                continue
            }
            marker[tmp_root] = 2;
            deg = 0;
            let mut search_size = 0;
            sn_size = 0;
            j_start = x_adj[tmp_root];
            j_stop = x_adj[tmp_root + 1];
            for j in j_start..j_stop {
                let mut node = adj[j];
                if marker[node] != 0 {
                    continue
                }
                if degree[node] >= 0 {
                    reach[reach_size + search_size] = node;
                    search_size += 1;
                    deg += q_size[node];
                    marker[node] = 1;
                    continue
                }
                neighbor[neighbor_size + sn_size] = node;
                sn_size += 1;
                marker[node] = -1;
                let mut k = x_adj[node];
                let mut k_stop = x_adj[node + 1];
                while k < k_stop {
                    node = adj[k];
                    k += 1;
                    if node >= num_vertex {
                        node -= num_vertex;
                        if node < num_vertex {
                            k = x_adj[node];
                            k_stop = x_adj[node + 1];
                            continue
                        }
                        break
                    }
                    if marker[node] != 0 {
                        continue
                    }
                    reach[reach_size + search_size] = node;
                    search_size += 1;
                    deg += q_size[node];
                    marker[node] = 1;
                }
            }

            degree[tmp_root] = set_power + deg - 1;
            counter += 1;
            if search_size != 0 {
                for j in 0..search_size {
                    marker[reach[reach_size + j]] = 0;
                }
            }

            if sn_size != 0 {
                for j in 0..sn_size {
                    marker[neighbor[neighbor_size + j]] = 0;
                }
            }
        }

        // Update thresh
        marker[root] = 0;
        for j in 0..reach_size {
            let index = reach[j];
            if marker[index] < 0 {
                continue
            }
            marker[index] = 0;
            let n_deg = degree[index];
            if n_deg < min_deg {
                min_deg = n_deg
            }
            if n_deg > thresh {
                continue
            }
            min_deg = thresh;
            thresh = n_deg;
            search = inv[index];
        }

        // Transform graph
        if neighbor_size != 0 {
            let mut nx_node = q_head[root];
            counter = 0;
            for i in 0..neighbor_size {
                let node = neighbor[i];
                q_link[nx_node] = node;
                counter += q_size[node];
                nx_node = q_head[node];
            }
            let mut node = root;
            q_head[root] = nx_node;
            q_size[root] += counter;
            j_start = x_adj[root];
            j_stop = x_adj[root + 1] - 1;
            for i in 0..reach_size {
                let index = reach[i];
                if marker[index] < 0 {
                    continue
                }
                if j_start == j_stop {
                    let l_node = j_start;
                    loop {
                        node = q_link[node];
                        j_start = x_adj[node];
                        j_stop = x_adj[node+1] - 1;
                        if j_start < j_stop {
                            break
                        }
                    }
                    adj[l_node] = node + num_vertex;
                }
                adj[j_start] = index;
                j_start += 1;
            }
            adj[j_start] = usize::MAX;
            for i in 0..reach_size {
                node = reach[i];
                if marker[node] < 0 {
                    continue
                }
                j_start = x_adj[node];
                j_stop = x_adj[node+1];
                for j in j_start..j_stop {
                    if marker[adj[j]] < 0 {
                        adj[j] = root;
                        break
                    }
                }
            }
        }
    }
}

fn merge(a: &mut Vec<usize>, a_len: usize, b: &mut [usize], b_len: usize, t: &mut Vec<usize>) -> usize {
    let mut i = 0_usize;
    let mut j = 0_usize;
    let mut k = 0_usize;

    while i < a_len && j < b_len {
        if a[i] < b[j] {
            t[k] = a[i];
            i += 1;
        } else if a[i] == b[j] {
            t[k] = a[i];
            i += 1;
            j += 1;
        } else {
            t[k] = b[j];
            j += 1;
        }
        k += 1;
    }

    while i < a_len {
        t[k] = a[i];
        k += 1;
        i += 1;
    }

    while j < b_len {
        t[k] = b[j];
        k += 1;
        j += 1;
    }
    k
}
fn sym_factorize(factor : &mut BccFactor, a_ptr: &Vec<usize>, a_ind:  &Vec<usize>) {
    let mut i_space = 0;
    factor.x_val.resize(factor.num_vertex + 1, 0);
    factor.x_ind.resize(factor.num_vertex, 0);
    
    // Allocate work memory
    let mut mrg_lnk = vec![0; factor.num_vertex];
    let mut stack = vec![0; factor.num_vertex];
    let mut node_set = vec![0; factor.num_vertex];
    
    // Allocate initial buffer for compressed row indices
    let mut l_size = 20 * factor.num_vertex;
    factor.l_ind.resize(l_size, 0);
    
    factor.x_val[0] = 0;
    
    // Set all lists as empty
    for i in 0..factor.num_vertex {
        mrg_lnk[i] = factor.num_vertex;
    }
    
    for k in 0..factor.num_vertex {
        let mut node = factor.perm[k];
        let start = a_ptr[node];
        let stop = a_ptr[node + 1];

        // Isolated root
        if start == stop {
            factor.x_ind[k] = i_space;
            factor.x_val[k + 1] = factor.l_space;
            continue
        }

        // Consolidate A(*,k)
        let mut set_size = 0;
        for j in start..stop {
            node = factor.inv_p[a_ind[j]];
            if node <= k {
                continue
            }
            node_set[set_size] = node;
            set_size += 1;
        }
        // sort.Ints(node_set[0:set_size])
        (& mut node_set[0..set_size]).sort();


        // Inspect a kids list and consolidate L(*,k)
        let mut i = mrg_lnk[k];
        while i < factor.num_vertex {
            // erge two sets
            set_size = merge(&mut node_set, set_size, &mut factor.l_ind[factor.x_ind[i] + 1..], factor.x_val[i + 1] - factor.x_val[i] - 1, &mut stack);
            // Now a set is in stack
            mem::swap(&mut stack, &mut node_set);
            i = mrg_lnk[i];
        }

        // Maybe compression?
        if k > 0 && factor.l_ind[factor.x_ind[k - 1]] == k && set_size == (factor.x_val[k] - factor.x_val[k - 1] - 1) {
            // It's an indistinguishable node
            factor.x_ind[k] = factor.x_ind[k - 1] + 1;
            factor.l_space += set_size;
            factor.x_val[k+1] = factor.l_space;
        } else {
            if i_space+set_size > l_size {
                l_size *= 2;
                factor.l_ind.resize(l_size, 0);
            }
            let mut j = i_space;
            for i in 0..set_size {
                factor.l_ind[j] = node_set[i];
                j += 1;
            }
            factor.x_ind[k] = i_space;
            i_space += set_size;
            factor.l_space += set_size;
            factor.x_val[k + 1] = factor.l_space;
        }
        if set_size > 1 {
            let i = factor.l_ind[factor.x_ind[k]];
            mrg_lnk[k] = mrg_lnk[i];
            mrg_lnk[i] = k;
        }
    }
}

fn sp_factor(factor : &mut BccFactor, matrix: &BccSparseMatrix, tolerance: f64) -> Result<(), FemError> {
    // Allocate main l_val array
    let d_enter = (matrix.block_size * (matrix.block_size + 1)) / 2;
    let mem_size = factor.l_space * matrix.block_size * matrix.block_size + matrix.num_vertex * d_enter;

    factor.d_val.resize(d_enter * matrix.num_vertex, 0.0);
    factor.l_val.resize(mem_size - d_enter * matrix.num_vertex, 0.0);

    factorize(factor, matrix, tolerance)
}

// ************************************************************
// GENERAL SPARSE SYMMETRIC SCHEME, CLASSIC FUN-IN ALGORITHM
// ************************************************************
fn factorize(factor : &mut BccFactor, matrix: &BccSparseMatrix, tolerance: f64) -> Result<(), FemError> {
    const MAX_BLOCK_SIZE: usize = 16;

    let mut first = vec![0; matrix.num_vertex];
    let mut link = vec![matrix.num_vertex; matrix.num_vertex];
    // Allocate diagonal flags array
    let mut s_val = vec![0; matrix.num_vertex * matrix.block_size];

    let mut translate = vec![0; matrix.num_vertex];
    let mut dd = vec![0.0; MAX_BLOCK_SIZE * MAX_BLOCK_SIZE];
    let mut hh = vec![0.0; MAX_BLOCK_SIZE * MAX_BLOCK_SIZE];
    let mut rd = vec![0.0; MAX_BLOCK_SIZE];

    // Sizes of diagonal & sub-diagonal blocks
    let f_enter = matrix.block_size * matrix.block_size;
    let d_enter = (matrix.block_size * (matrix.block_size + 1)) / 2;

    // Setup initial diagonal range
    let mut d_range: f64;
    let mut knew: usize;

    let mut msg = Messenger::new("Factorization equations", 0, factor.num_vertex as i64, 1);
    // Main routine
    for j in 0..matrix.num_vertex {
        // if *canceled {
        //     return fmt.Errorf("the process was stopped by user")
        // }
        msg.add_progress();
        let mut start = factor.x_val[j];
        let mut stop = factor.x_val[j + 1];
        // Set entries into a factored column
        let mut i_sub = factor.x_ind[j];
        let mut offset = f_enter * start;
        for _ in start..stop {
            translate[factor.l_ind[i_sub]] = offset;
            i_sub += 1;
            for k in 0..f_enter {
                factor.l_val[offset + k] = 0.0;
            }
            offset += f_enter;
        }

        // Form A structure
        let mut node = factor.perm[j];
        start = matrix.ptr[node];
        stop = matrix.ptr[node + 1];

        offset = f_enter * start;
        for i in start..stop {
            node = factor.inv_p[matrix.ind[i]];
            if node >= j {
                for k in 0..f_enter {
                    if node == j {
                        dd[k] = matrix.val[offset + k];
                    } else {
                        factor.l_val[translate[node]+k] = matrix.val[offset + k];
                    }
                }
            }
            offset += f_enter
        }

        // Merge with previous columns
        let mut k = link[j];
        while k < matrix.num_vertex {
            knew = link[k];
            start = first[k];
            // Calculate HH
            offset = 0;
            for bi in 0..matrix.block_size {
                let value = s_val[matrix.block_size*k+bi] as f64;
                for bj in 0..matrix.block_size {
                    hh[offset+bj] = factor.l_val[f_enter * start + offset + bj] * value;
                }
                offset += matrix.block_size;
            }

            // DD -= LL * HH
            offset = 0;
            for bi in 0..matrix.block_size {
                for bj in 0..matrix.block_size {
                    let value = hh[matrix.block_size * bj + bi];
                    for bk in bi..matrix.block_size {
                        dd[offset + bk] -= factor.l_val[f_enter * start + matrix.block_size * bj + bk] * value;
                    }
                }
                offset += matrix.block_size;
            }
            start += 1;
            stop = factor.x_val[k + 1];
            if start < stop {
                // Update links
                first[k] = start;
                let mut i = factor.x_ind[k] + start - factor.x_val[k];
                i_sub = factor.l_ind[i];
                link[k] = link[i_sub];
                link[i_sub] = k;

                // Merge sub-diagonal
                for n in start..stop {
                    for bi in 0..matrix.block_size {
                        for bj in 0..matrix.block_size {
                            let value = hh[matrix.block_size * bj + bi];
                            for bk in 0..matrix.block_size {
                                factor.l_val[translate[factor.l_ind[i]]+bk+matrix.block_size*bi] -= factor.l_val[f_enter*start+matrix.block_size*bj+bk+f_enter*(n-start)] * value;
                            }
                        }
                    }
                    i += 1;
                }
            }

            k = knew;
        }

        // Diagonal block factorization
        offset = 0;
        for bi in 0..matrix.block_size {
            // Accumulating
            for bj in 0..bi {
                let value = dd[matrix.block_size*bj+bi] * s_val[matrix.block_size*j+bj] as f64;
                for bk in bi..matrix.block_size {
                    dd[offset + bk] -= dd[matrix.block_size * bj + bk] * value;
                }
            }
            let mut value = dd[offset + bi];

            // Check indefinite
            if value < 0.0 {
                return Err(FemError::SingularMatrixError);
            } else {
                d_range = 1.0;
                s_val[matrix.block_size * j + bi] = 1;
            }

            // Check numerical stability
            if value < tolerance {
                return Err(FemError::SingularMatrixError);
            }

            value = value.sqrt();
            dd[offset + bi] = value;
            value = d_range / value;
            rd[bi] = value;

            // Scale sub-diagonal
            for bk in bi + 1..matrix.block_size {
                dd[offset + bk] *= value;
            }
            offset += matrix.block_size;
        }

        // Save factored diagonal block
        let mut col_size = matrix.block_size;
        offset = d_enter * j;
        let mut shift = 0;
        for bi in 0..matrix.block_size {
            for bj in bi..matrix.block_size {
                factor.d_val[offset + bj] = dd[shift + bj];
            }
            col_size -= 1;
            offset += col_size;
            shift += matrix.block_size;
        }

        // Do we must factor sub-diagonal?
        start = factor.x_val[j];
        stop = factor.x_val[j + 1];
        if start >= stop {
            continue
        }

        // Update links
        first[j] = start;
        let i = factor.x_ind[j];
        i_sub = factor.l_ind[i];
        link[j] = link[i_sub];
        link[i_sub] = j;

        // Calculate HH
        offset = 0;
        shift = 0;
        for bi in 0..matrix.block_size {
            let value = s_val[matrix.block_size * j + bi] as f64;
            for bj in bi + 1..matrix.block_size {
                hh[shift + bj] = dd[offset + bj] * value;
            }
            offset += matrix.block_size;
            shift += matrix.block_size;
        }

        // Scale sub-diagonal
        offset = f_enter * start;
        for _ in start..stop {
            shift = offset;
            for bi in 0..matrix.block_size {
                for bj in 0..bi {
                    let value = hh[matrix.block_size * bj + bi];
                    for bk in 0..matrix.block_size {
                        factor.l_val[shift + bk] -= factor.l_val[offset + matrix.block_size * bj + bk] * value;
                    }
                }
                let value = rd[bi];
                for bk in 0..matrix.block_size {
                    factor.l_val[shift + bk] *= value;
                }
                shift += matrix.block_size;
            }
            offset += f_enter;
        }
    }

    Ok(())
}

fn sp_solve(factor : &mut BccFactor, vector: &mut Vec<f64>) {
    let mut msg = Messenger::new("Solution of the system of equations", 0, 0, 0);
    // Direct permutation
    permutation(vector, &factor.perm, factor.num_vertex, factor.block_size);
    // Choose an optimized variant or universal
    sparse_solve(factor, vector);
    // Invert permutation
    permutation(vector, &factor.inv_p, factor.num_vertex, factor.block_size);
    msg.stop();
}

fn permutation(rhs: &mut Vec<f64>, order: &Vec<usize>, num_vertex: usize, block_size: usize) {
    for i in 0..num_vertex {
        let mut node = order[i];
        if node == i {
            continue
        }
        while node < i {
            node = order[node];
        }
        for j in 0..block_size {
            (rhs[block_size * node + j], rhs[block_size * i + j]) = (rhs[block_size * i + j], rhs[block_size * node + j])
        }
    }
}

// ******************************************************************
//
//	GENERAL SPARSE SYMMETRIC SYSTEM
//
// ******************************************************************
fn sparse_solve(factor : &mut BccFactor, rhs: &mut Vec<f64>) {
    let f_enter = factor.block_size * factor.block_size;
    let d_enter = (factor.block_size * (factor.block_size + 1)) / 2;
    
    // Forward substitution
    let mut offset = 0;
    for j in 0..factor.num_vertex {
    
        // Dense lover diagonal system
        let mut col_ize = factor.block_size;
        let mut run = 0;
        let mut shift = 0;
        for col in 0..factor.block_size {
            let mut value = rhs[col + j * factor.block_size];
            if value != 0.0 {
                value /= factor.d_val[offset + shift + col];
                rhs[col + j * factor.block_size] = value;
                run = 1;
                for row in col + 1..factor.block_size {
                    rhs[row+j*factor.block_size] -= factor.d_val[offset + shift + row] * value;
                }
            }
            col_ize -= 1;
            shift += col_ize;
        }

        if run != 0 {
            let start = factor.x_val[j];
            let stop = factor.x_val[j + 1];
            if start < stop {
                let mut i = factor.x_ind[j];
                shift = f_enter * start;
                for _ in start..stop {
                    for col in 0..factor.block_size {
                        let value = rhs[col + j * factor.block_size];
                        for row in 0..factor.block_size {
                            rhs[row + factor.block_size * factor.l_ind[i]] -= factor.l_val[shift + factor.block_size * col + row] * value;
                        }
                    }
                    i += 1;
                    shift += f_enter;
                }
            }
        }
        offset += d_enter;
    }
    
    // Backward substitution
    let mut offset:i32 = (factor.num_vertex * d_enter - factor.block_size) as i32;
    for j in (0..factor.num_vertex).rev() {
        let start = factor.x_val[j];
        let stop = factor.x_val[j + 1];

        if start < stop {
            let mut i = factor.x_ind[j];
            let mut shift = f_enter * start;
            for _ in start..stop {
                for row in 0..factor.block_size {
                    let mut value = 0.0;
                    for col in 0..factor.block_size {
                        value += factor.l_val[shift + row * factor.block_size + col] * rhs[col + factor.block_size * factor.l_ind[i]];
                    }
                    rhs[row+j*factor.block_size] -= value;
                }
                i += 1;
                shift += f_enter;
            }
        }

        // Dense upper diagonal system
        let mut col_size = 0;
        let mut shift = 0i32;
        for row in (0..factor.block_size).rev() {
            let mut value = 0.0;
            for col in row + 1..factor.block_size {
                let index :i32 = col as i32 + offset + shift;
                value += rhs[col + j * factor.block_size] * factor.d_val[index as usize];
            }
            let index :i32 = row as i32 + offset + shift;
            rhs[row + j * factor.block_size] = (rhs[row + j * factor.block_size] - value) / factor.d_val[index as usize];
            col_size += 1;
            shift -= col_size;
        }
        offset -= d_enter as i32;
    }
}
