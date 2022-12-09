#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <deque>
#include <unordered_set>
#include <atomic>
#include <chrono>

#include "colorSCC.hpp"
#include "sparse_util.hpp"

#include <omp.h>

#define UNCOMPLETED_SCC_ID 18446744073709551615
#define MAX_COLOR 18446744073709551615

// For the first time only, where all SCC_ids are -1
size_t trimVertices_inplace_first_time(const Sparse_matrix& inb, const Sparse_matrix& onb, std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    //std::atomic<size_t> trimed(0);
    std::atomic<size_t> trimed(0);

    # pragma omp parallel for shared(trimed)
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

// first time only, but only one matrix
size_t trimVertices_inplace_first_time_single_direction(const Sparse_matrix& nb, std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    std::atomic<size_t> trimed(0);
    //std::vector<bool> hasOtherWay(nb.n, false);

    std::deque<std::atomic<bool>> hasOtherWay;
    for(size_t i = 0; i < nb.n; i++) {
        hasOtherWay.emplace_back(false);
    }

 //   bool hasOtherWay[nb.n];
    //std::fill(hasOtherWay, hasOtherWay + nb.n, false);

    # pragma omp parallel for shared(trimed, hasOtherWay)
    for(size_t source = 0; source < nb.n; source++) {
        if(nb.ptr[source] == nb.ptr[source + 1]) {
            SCC_id[source] = SCC_count + ++trimed;
        } else {
            for (size_t i = nb.ptr[source]; i < nb.ptr[source + 1]; i++) {
                size_t neighbor = nb.val[i];
                hasOtherWay[neighbor] = true;
            }
        }
    } 

    # pragma omp parallel for shared(trimed)
    for(size_t source = 0; source < nb.n; source++) {
        if(!hasOtherWay[source] && SCC_id[source] == UNCOMPLETED_SCC_ID) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }

    return trimed;
}

// vleft + onb
size_t trimVertices_inplace(const Sparse_matrix& inb, const Sparse_matrix& onb, const std::vector<size_t>& vleft,
                                    std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    std::atomic<size_t> trimed(0);
    //size_t trimed = 0;

    # pragma omp parallel for shared(trimed)
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

    return trimed;
}

size_t trimVertices_inplace_single_direction(const Sparse_matrix& nb, const std::vector<size_t>& vleft,
                                        std::vector<size_t>& SCC_id, size_t SCC_count) { 

    //size_t trimed = 0;
    std::atomic<size_t> trimed(0);
    const size_t vertices_left = vleft.size();
    const size_t n = nb.n;

    std::deque<std::atomic<bool>> hasOtherWay;
    for(size_t i = 0; i < n; i++) {
        hasOtherWay.emplace_back(false);
    }

    // only going on the vertices left, no need to check for scc_id
    #pragma omp parallel for shared(trimed, hasOtherWay)
    for(size_t index = 0; index < vertices_left; index++) {
        size_t source = vleft[index];

        bool hasOneWay = false;
        for(size_t i = nb.ptr[source]; i < nb.ptr[source + 1]; i++) {
            size_t neighbor = nb.val[i];

            // if SCC_id[neighbor] == UNCOMPLETED_SCC_ID, then neighbor in vleft
            if(SCC_id[neighbor] == UNCOMPLETED_SCC_ID) {
                hasOneWay = true;
                // need index of neighbor
                hasOtherWay[neighbor] = true;
            }
        }

    // no inc neighbors then surely trim
        if(!hasOneWay) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }

    #pragma omp parallel for shared(trimed)
    for(size_t source = 0; source < n; source++) {
        // hasOtherWayIndex] == false => noone in vleft points/is pointed to (depending on if nb is incoming or outgoing) to source
        if(!hasOtherWay[source] && SCC_id[source] == UNCOMPLETED_SCC_ID) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    return trimed;
}

void bfs_colors_inplace( const Sparse_matrix& nb, const size_t source, std::vector<size_t>& SCC_id,
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

std::vector<size_t> colorSCC(Coo_matrix& M, bool DEBUG) {
    Sparse_matrix inb;
    Sparse_matrix onb;

    DEB("Starting conversion");
    
    # pragma omp parallel sections
    {
        # pragma omp section
        {
            coo_tocsr(M, inb);
        }
        # pragma omp section
        {
            coo_tocsc(M, onb);
        }
    }
    // if we are poor on memory, we can free M
    M.Ai = std::vector<size_t>();
    M.Aj = std::vector<size_t>();
    DEB("Finished conversion");

    return colorSCC_no_conversion(inb, onb, true, DEBUG);
}

// If !USE_ONB then we never use onb, but we still need to pass it, an optional wouldn't work because we need to pass it by reference
std::vector<size_t> colorSCC_no_conversion(const Sparse_matrix& inb, const Sparse_matrix& onb, bool USE_ONB, bool DEBUG) {
    size_t n = inb.n;

    std::vector<size_t> SCC_id(n, UNCOMPLETED_SCC_ID);
    size_t SCC_count = 0;

    std::vector<size_t> vleft(n);
    for (size_t i = 0; i < n; i++) {
        vleft[i] = i;
    }

    DEB("First time trim")
    if(USE_ONB) {
        SCC_count += trimVertices_inplace_first_time(inb, onb, SCC_id, SCC_count);
    } else {
        SCC_count += trimVertices_inplace_first_time_single_direction(inb, SCC_id, SCC_count);
    }
    DEB("Finished trim")

    DEB("First erasure")
    std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
    DEB("Finished first erasure")
    DEB("Size difference: " << SCC_count)

    std::vector<size_t> colors(n);

    size_t iter = 0;
    size_t total_tries = 0;
    while(!vleft.empty()) {
        iter++;
        DEB("Starting while loop iteration " << iter)

        # pragma omp parallel for shared(colors)
        for(size_t i = 0; i < n; i++) {
            colors[i] = SCC_id[i] == UNCOMPLETED_SCC_ID ? i : MAX_COLOR;
        }

        DEB("Starting to color")
        bool made_change = true;
        while(made_change) {
            made_change = false;

            total_tries++;
            # pragma omp parallel for shared(colors, made_change, vleft)
            for(size_t i = 0; i < vleft.size(); i++) {
                size_t u = vleft[i];

                for(size_t j = inb.ptr[u]; j < inb.ptr[u + 1]; j++) {
                    size_t v = inb.val[j];

                    size_t new_color = colors[v];

                    // if the neightbor v is in some SCC, then it's color will be MAX_COLOR hance not triggering the if
                    if(new_color < colors[u]) {
                        colors[u] = new_color;

                        made_change = true;
                    }
                }
            }
            iter += vleft.size();
        }
        DEB("Finished coloring")

        DEB("Set of colors part")
        auto unique_colors_set = std::unordered_set<size_t> (colors.begin(), colors.end());
        unique_colors_set.erase(MAX_COLOR);

        DEB("Found " << unique_colors_set.size() << " unique colors")

        auto unique_colors = std::vector<size_t>(unique_colors_set.begin(), unique_colors_set.end());
        DEB("End of Set of colors part");

        DEB("Starting bfs")
        # pragma omp parallel for 
        for(size_t i = 0; i < unique_colors.size(); i++) {
            const size_t color = unique_colors[i];
            const size_t _SCC_count = SCC_count + i + 1;

            bfs_colors_inplace(inb, color, SCC_id, _SCC_count, colors, color);
        }
        SCC_count += unique_colors.size();
        DEB("Finished BFS")

        DEB("Trim + erasure")
        // remove all vertices that are in some SCC
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });

        if (USE_ONB) {
           SCC_count += trimVertices_inplace(inb, onb, vleft, SCC_id, SCC_count);
        } else {
           SCC_count += trimVertices_inplace_single_direction(inb, vleft, SCC_id, SCC_count);
        }
        // clean up vleft after trim
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
        DEB("Finished trim + erasure")

    }
    DEB("Total tries: " << total_tries)
    DEB("Total iterations: " << iter)
    DEB("Total SCCs: " << SCC_count)
    return SCC_id;
}
