#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <set>
#include <unordered_set>
#include <chrono>

#include "sparse_util.hpp"
#include "colorSCC.hpp"

#include <coz.h>

#define DEB(x) if(DEBUG) {std::cout << x << std::endl;}

#define UNCOMPLETED_SCC_ID 18446744073709551615
#define MAX_COLOR 18446744073709551615

// For the first time only, where all SCC_ids are -1
size_t trimVertices_inplace_normal_first_time(const Sparse_matrix& inb, const Sparse_matrix& onb, std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    size_t trimed = 0;

    for(size_t source = 0; source < inb.n; source++) {
        bool hasIncoming = inb.ptr[source] != inb.ptr[source + 1];

        bool hasOutgoing = onb.ptr[source] != onb.ptr[source + 1];

        if(!hasIncoming | !hasOutgoing) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    //std::cout << "trimed: " << trimed << std::endl;

    return trimed;
}

// vleft + onb
size_t trimVertices_inplace_normal(const Sparse_matrix& inb, const Sparse_matrix& onb, const std::vector<size_t>& vleft,
                                    std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    size_t trimed = 0;

    for(size_t index = 0; index < vleft.size(); index++) {
        const size_t source = vleft[index];

        bool hasIncoming = false;
        for(size_t i = inb.ptr[source]; i < inb.ptr[source + 1]; i++) {
            if(SCC_id[inb.val[i]] == UNCOMPLETED_SCC_ID) {
                hasIncoming = true;
                break;
            }
        }

        bool hasOutgoing = false;
        for(size_t i = onb.ptr[source]; i < onb.ptr[source + 1]; i++) {
            if(SCC_id[onb.val[i]] == UNCOMPLETED_SCC_ID) {
                hasOutgoing = true;
                break;
            }
        }

        if(!hasIncoming | !hasOutgoing) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    //std::cout << "trimed: " << trimed << std::endl;

    return trimed;
}

size_t trimVertices_inplace_normal_no_onb(const Sparse_matrix& inb, const std::vector<size_t>& vleft,
                                        std::vector<size_t>& SCC_id, size_t SCC_count) { 

    size_t trimed = 0;
    const size_t vertices_left = vleft.size();
    const size_t n = inb.n;

    auto hasOutgoing = std::vector<bool>(n, false);

    for(size_t index = 0; index < vertices_left; index++) {
        size_t source = vleft[index];

        bool hasIncoming = false;
        for(size_t i = inb.ptr[source]; i < inb.ptr[source + 1]; i++) {
            size_t neighbor = inb.val[i];

            // if SCC_id[neighbor] == UNCOMPLETED_SCC_ID, then neighbor in vleft
            if(SCC_id[neighbor] == UNCOMPLETED_SCC_ID) {
                hasIncoming = true;
                // need index of neighbor
                hasOutgoing[neighbor] = true;
            }
        }

    // no inc neighbors then surely trim
        if(!hasIncoming) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }

    for(size_t source = 0; source < n; source++) {
        // check if it has already been trimmed in the prev step
        if (SCC_id[source] != UNCOMPLETED_SCC_ID) continue;

        // noone in vleft was pointed to by source, so source is surely trimable
        // hasOutgoing[index] == false => noone in vleft points to source = vleft[index]
        if(!hasOutgoing[source]) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    return trimed;
}

void bfs_sparse_colors_all_inplace( const Sparse_matrix& nb, const size_t source, std::vector<size_t>& SCC_id,
                                const size_t SCC_count, const std::vector<size_t>& colors, const size_t color) {
    SCC_id[source] = SCC_count;

    std::queue<size_t> q;
    q.push(source);

    while(!q.empty()) {
        size_t v = q.front();
        q.pop();

        for(size_t i = nb.ptr[v]; i < nb.ptr[v + 1]; i++) {
            size_t u = nb.val[i];

            if(SCC_id[u] == UNCOMPLETED_SCC_ID && colors[u] == color) {
                SCC_id[u] = SCC_count;
                q.push(u);
            }
        }
    }
}

std::vector<size_t> colorSCC_no_conversion(const Sparse_matrix& inb, const Sparse_matrix& onb, bool DEBUG);

std::vector<size_t> colorSCC(Coo_matrix& M, bool DEBUG) {
    Sparse_matrix inb;
    Sparse_matrix onb;

    DEB("Starting conversion");
    
    COZ_BEGIN("convert");

    coo_tocsr(M, inb);
    coo_tocsc(M, onb);

    COZ_END("convert");
    // if we are poor on memory, we can free M
    M.Ai = std::vector<size_t>();
    M.Aj = std::vector<size_t>();
    DEB("Finished conversion");

    return colorSCC_no_conversion(inb, onb, DEBUG);
}

// working
std::vector<size_t> colorSCC_no_conversion(const Sparse_matrix& inb, const Sparse_matrix& onb, bool DEBUG) {
    size_t n = inb.n;

    std::vector<size_t> SCC_id(n);
    std::fill(SCC_id.begin(), SCC_id.end(), UNCOMPLETED_SCC_ID);
    size_t SCC_count = 0;


    std::vector<size_t> vleft(n);
    for (size_t i = 0; i < n; i++) {
        vleft[i] = i;
    }

    DEB("Starting trim")
    SCC_count += trimVertices_inplace_normal_first_time(inb, onb, SCC_id, SCC_count);
    DEB("Finished trim")

    std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });

    DEB("Size difference: " << SCC_count)

    std::vector<size_t> colors(n);

    size_t iter = 0;
    size_t total_tries = 0;
    while(!vleft.empty()) {
        iter++;

        DEB("Starting while loop iteration " << iter)

        for(size_t i = 0; i < n; i++) {
            colors[i] = (SCC_id[i] == UNCOMPLETED_SCC_ID) ? i : MAX_COLOR;
        }

        DEB("Starting to color")
        bool made_change = true;
        while(made_change) {
            made_change = false;

            total_tries++;
            for(size_t i = 0; i < vleft.size(); i++) {
                size_t u = vleft[i];

                for(size_t i = inb.ptr[u]; i < inb.ptr[u + 1]; i++) {
                    size_t v = inb.val[i];
                    //if(colors[v] == MAX_COLOR) continue;
                    size_t new_color = colors[v];

                    if(new_color < colors[u]) {
                        colors[u] = new_color;
                        made_change = true;
                    }
                }
            }
        }
        DEB("Finished coloring")

        auto unique_colors_set = std::unordered_set<size_t> (colors.begin(), colors.end());
        unique_colors_set.erase(MAX_COLOR);

        DEB("Found " << unique_colors_set.size() << " unique colors")

        auto unique_colors = std::vector<size_t>(unique_colors_set.begin(), unique_colors_set.end());

        //std::sort(unique_colors.begin(), unique_colors.end());

        DEB("Starting BFS")
        for(size_t i = 0; i < unique_colors.size(); i++) {
            size_t color = unique_colors[i];
            const size_t _SCC_count = SCC_count + i + 1;

            bfs_sparse_colors_all_inplace(inb, color, SCC_id, _SCC_count, colors, color);
        }
        DEB("Finished BFS")
        SCC_count += unique_colors.size();

        // remove all vertices that are in some SCC
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
        SCC_count += trimVertices_inplace_normal(inb, onb, vleft, SCC_id, SCC_count);

        // clean up vleft after trim
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
 
    }

    std::cout << "Total tries: " << total_tries << std::endl;
    std::cout << "Total iterations: " << iter << std::endl;


    return SCC_id;
}
