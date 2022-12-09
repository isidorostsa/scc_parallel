#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <queue>
#include <unordered_set>
#include <chrono>
#include <iterator>
#include <list>
#include <optional>
#include <atomic>

#include "sparse_util.hpp"
#include "colorSCC.hpp"

#include <pthread.h>

#define BFS_GRAIN_SIZE 3
#define COLOR_GRAIN_SIZE 1000

#define UNCOMPLETED_SCC_ID -1
#define MAX_COLOR -1

/**
 * @brief A trimming function, without checking for removed vertices, and without taking into account the vertices it trims in the neighbor number calculation
 * @param inb incoming neighbors
 * @param onb outgoing neighbors
 * @param SCC_id the SCC id of each vertex, needed to save assign the SCC id of the trimmed vertices
 * @param SCC_count the number of SCCs found so far, needed to calculate the scc id of the trimmed vertices
 * @return the number of trimmed vertices
 */size_t trimVertices_inplace_first_time(const Sparse_matrix& inb, const Sparse_matrix& onb, std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    std::atomic<size_t> trimed(0);

    for(size_t source = 0; source < inb.n; source++) {
        // the distance between the pointers is the number of neighbors
        // 0 means no neighbors
        bool hasIncoming = inb.ptr[source] != inb.ptr[source + 1];

        bool hasOutgoing = onb.ptr[source] != onb.ptr[source + 1];

        if(!hasIncoming | !hasOutgoing) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    return trimed;
}

/**
 * @brief A trimming function, without checking for removed vertices, and with only knowing the neighbors in one direction
 * @param nb neighbors
 * @param SCC_id the SCC id of each vertex, needed to save assign the SCC id of the trimmed vertices
 * @param SCC_count the number of SCCs found so far, needed to calculate the scc id of the trimmed vertices
 * @return the number of trimmed vertices
 */
size_t trimVertices_inplace_first_time_single_direction(const Sparse_matrix& nb, std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    size_t trimed = 0;

    std::vector<bool> hasOtherWay(nb.n, false);

    for(size_t source = 0; source < nb.n; source++) {
        // the distance between the pointers is the number of neighbors
        if(nb.ptr[source] == nb.ptr[source + 1]) {
            SCC_id[source] = SCC_count + ++trimed;
        } else {
            for (size_t i = nb.ptr[source]; i < nb.ptr[source + 1]; i++) {
                size_t neighbor = nb.val[i];
                hasOtherWay[neighbor] = true;
            }
        }
    } 

    // trim the vertices are not neighbors of any other vertex in the other direction, and have not been trimmed yet
    for(size_t source = 0; source < nb.n; source++) {
        if(!hasOtherWay[source] && SCC_id[source] == UNCOMPLETED_SCC_ID) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }

    return trimed;
}

/**
 * @brief A trimming function. Assumes that the vertices in vleft are not trimmed yet, and that the SCC_id of the trimmed vertices is already set. Changes SCC_IDs, does not change vleft.
 * @param inb incoming neighbors
 * @param onb outgoing neighbors 
 * @param vleft the vertices that are not trimmed yet
 * @param SCC_id the SCC id of each vertex
 * @param SCC_count the number of SCCs found so far, needed to calculate the scc id of the trimmed vertices
 * @return the number of trimmed vertices
 */
size_t trimVertices_inplace(const Sparse_matrix& inb, const Sparse_matrix& onb, const std::vector<size_t>& vleft,
                                    std::vector<size_t>& SCC_id, const size_t SCC_count) { 
    size_t trimed = 0;

    // check for non-trimmed neighbors of the source vertex in both directions
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

/**
 * @brief A trimming function. Assumes that the vertices in vleft are not trimmed yet, and that the SCC_id of the trimmed vertices is already set. Changes SCC_IDs, does not change vleft.
 * @param nb neighbors
 * @param vleft the vertices that are not trimmed yet
 * @param SCC_id the SCC id of each vertex
 * @param SCC_count the number of SCCs found so far, needed to calculate the scc id of the trimmed vertices
 * @return the number of trimmed vertices
 */
size_t trimVertices_inplace_single_direction(const Sparse_matrix& nb, const std::vector<size_t>& vleft,
                                        std::vector<size_t>& SCC_id, size_t SCC_count) { 

    //size_t trimed = 0;
    size_t trimed = 0;
    const size_t vertices_left = vleft.size();
    const size_t n = nb.n;

    auto hasOtherWay = std::vector<bool>(n, false);

    // only going on the vertices left, no need to check for scc_id
    for(size_t index = 0; index < vertices_left; index++) {
        size_t source = vleft[index];

        bool hasOneWay = false;
        for(size_t i = nb.ptr[source]; i < nb.ptr[source + 1]; i++) {
            size_t neighbor = nb.val[i];

            // if SCC_id[neighbor] == UNCOMPLETED_SCC_ID, then neighbor in vleft
            if(SCC_id[neighbor] == UNCOMPLETED_SCC_ID) {
                hasOneWay = true;
                hasOtherWay[neighbor] = true;
            }
        }

    // no inc neighbors then surely trim
        if(!hasOneWay) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }

    for(size_t source = 0; source < n; source++) {
        // hasOtherWayIndex] == false => noone in vleft points/is pointed to (depending on if nb is incoming or outgoing) to source
        if(!hasOtherWay[source] && SCC_id[source] == UNCOMPLETED_SCC_ID) {
            SCC_id[source] = SCC_count + ++trimed;
        }
    }
    return trimed;
}


/**
 * @brief BFS that changes the SCC_id of all vertices it reaches that have the given color. Assumes that the SCC_id of the trimmed vertices is already set.
 * @param nb neighbors in the direction of the BFS
 * @param source the source vertex
 * @param SCC_id the SCC id of each vertex
 * @param SCC_count the number of SCCs found so far, needed to calculate the scc id of the assigned vertices
 * @param colors the color of each vertex
 * @param color the color of the vertices that will be assigned the SCC_count
 * @return (void)
 */
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

            if(colors[u] == color && SCC_id[u] != SCC_count) {
                SCC_id[u] = SCC_count;
                q.push(u);
            }
        }
    }
}

/**
 * @brief Needed for the parallel BFS. Contains all the information needed for the BFS.
 */
struct bfs_partitions_runner_struct
{
    const Sparse_matrix* nb;
    std::vector<size_t>* SCC_id;
    size_t SCC_count;

    // source and color are found from those below
    std::vector<size_t>* colors;
    std::vector<size_t>* unique_colors;
    size_t start;
    size_t end;

    bool should_quit_after;
};

/**
 * @brief Runs the bfs starting from a given set of sources. It represents one thread in the parallel BFS.
 * @param bfs_plus_info A struct containing all the information needed for the BFS
 * @return (void)
 */
void bfs_partitions_runner(bfs_partitions_runner_struct* bfs_plus_info) {
    const size_t SCC_count = bfs_plus_info->SCC_count;
    const size_t start = bfs_plus_info->start;
    const size_t end = bfs_plus_info->end;
    const Sparse_matrix& nb = *bfs_plus_info->nb;
    std::vector<size_t>& SCC_id = *bfs_plus_info->SCC_id;
    std::vector<size_t>& colors = *bfs_plus_info->colors;
    std::vector<size_t>& unique_colors = *bfs_plus_info->unique_colors;

    for(size_t i = start; i < end; i++) {
        size_t color = unique_colors[i];
        size_t _SCC_count = SCC_count + i + 1;
        bfs_colors_inplace(nb, color, SCC_id, _SCC_count, colors, color);
    }

    if(bfs_plus_info->should_quit_after) pthread_exit(0);
    return;
}

/**
 * @brief Needed for the parallel coloring. Contains all the information needed for the BFS.
 */
struct coloring_partitions_runner_struct
{
    const Sparse_matrix* inb;
    std::vector<size_t>* colors;
    const std::vector<size_t>* vleft;
    size_t start;
    size_t end;
    bool* made_change;
    bool should_quit_after;
};

/**
 * @brief Runs the coloring starting for a given set of sources. It represents one thread in the parallel coloring.
 * @param coloring_info A struct containing all the information needed for the coloring
 * @return (void)
 */
void coloring_partitions_runner(coloring_partitions_runner_struct* coloring_info) {
    const size_t start = coloring_info->start;
    const size_t end = coloring_info->end;

    const Sparse_matrix& inb = *coloring_info->inb;
    std::vector<size_t>& colors = *coloring_info->colors;
    const std::vector<size_t>& vleft = *coloring_info->vleft;

    bool& made_change = *coloring_info->made_change;

    for(size_t i = start; i < end; i++) {
        size_t u = vleft[i];

        for(size_t j = inb.ptr[u]; j < inb.ptr[u + 1]; j++) {
            size_t v = inb.val[j];

            size_t new_color = colors[v];
            if(new_color < colors[u]) {
                colors[u] = new_color;
                made_change = true;
            }
        }
    }
}

/**
 * @brief Finds the SCCs of a directed graph. Assumes that the graph is in csr and csc format. onb is optional.
 * @param inb incoming neighbors
 * @param onb outgoing neighbors (optional)
 * @param USE_ONB if true, onb is used, otherwise it is ignored
 * @param DEBUG if true, prints debug info
 * @return the SCC id of each vertex
 */
std::vector<size_t> colorSCC_no_conversion(const Sparse_matrix& inb, const Sparse_matrix& onb, bool USE_ONB, bool DEBUG, const size_t NUM_THREADS) {

    size_t n = inb.n;
    std::vector<size_t> SCC_id(n, UNCOMPLETED_SCC_ID);
    size_t SCC_count = 0;

    // an array of vertices that are left to be processed
    DEB("Starting trim")
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

    // remove all vertices that have been trimmed
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

        // the removed vertices are colored MAX_COLOR, this color never gets propagated to other vertices
        for(size_t i = 0; i < n; i++) {
            colors[i] = SCC_id[i] == UNCOMPLETED_SCC_ID ? i : MAX_COLOR;
        }

        DEB("Starting to color")
        bool made_change = true;
        while(made_change) {
            made_change = false;

            const size_t num_vertices = vleft.size();

            // Each thread gets at least COLOR_GRAIN_SIZE vertices to color
            // we may not use all threads
            size_t threads_to_use = (num_vertices/COLOR_GRAIN_SIZE <= 1) ? 1 : 
                                    (num_vertices/COLOR_GRAIN_SIZE > NUM_THREADS) ? NUM_THREADS : 1;

            DEB("Using " << threads_to_use << " threads for coloring.")

            std::vector<pthread_t> threads(threads_to_use);
            std::vector<coloring_partitions_runner_struct> coloring_info(threads_to_use);

            for(size_t i = 0; i < threads_to_use; i++) {

                // equally split partitions
                const size_t start = i * num_vertices / threads_to_use;
                const size_t end = (i == threads_to_use - 1) ? num_vertices : (i + 1) * num_vertices / threads_to_use;

                // info for each thread
                coloring_info[i] = {&inb, &colors, &vleft, start, end, &made_change, false};
                
                // a thread runs the coloring for all vertices in vleft[start:end]
                pthread_create(&threads[i], NULL, (void*(*)(void*))coloring_partitions_runner, &coloring_info[i]);
            }

            for(size_t i = 0; i < threads_to_use; i++) {
                pthread_join(threads[i], NULL);
            } 
        }

        DEB("Finished coloring")

        // Create a set of the unique colors to schedule the BFS
        // A BFS starts from each unique color
        DEB("Set of colors part");
        auto unique_colors_set = std::unordered_set<size_t> (colors.begin(), colors.end());
        unique_colors_set.erase(MAX_COLOR);
        DEB("Found " << unique_colors_set.size() << " unique colors")
        auto unique_colors = std::vector<size_t>(unique_colors_set.begin(), unique_colors_set.end());
        DEB("Set of colors part");

        DEB("Starting bfs")

        const size_t num_colors = unique_colors.size();
        // each thread will get a different set of colors to work on

        // Each thread gets at least BFS_GRAIN_SIZE vertices to color
        // we may not use all threads
        size_t threads_to_use = (num_colors/BFS_GRAIN_SIZE <= 1) ? 1 : 
                                (num_colors/BFS_GRAIN_SIZE > NUM_THREADS) ? NUM_THREADS : 1;

        DEB("Using " << threads_to_use << " threads for BFS")

        if(threads_to_use == 1) {
            bfs_partitions_runner_struct bfs_plus_info = 
                        {&inb, &SCC_id, SCC_count, &colors, &unique_colors, 0, num_colors, false};

            // no need to dispatch any threads
            bfs_partitions_runner(&bfs_plus_info);
        } else {

            std::vector<pthread_t> threads(threads_to_use);
            std::vector<bfs_partitions_runner_struct> bfs_plus_infos(threads_to_use);

            for(size_t i = 0; i < threads_to_use; i++) {
                // equally split partitions
                size_t start = i * num_colors / threads_to_use;
                size_t end = (i == threads_to_use - 1) ? num_colors : (i + 1) * num_colors / threads_to_use;

                // info for each thread
                // A thread starts a BFS from each color in unique_colors[start:end]
                bfs_plus_infos[i] = {&inb, &SCC_id, SCC_count, &colors, &unique_colors, start, end, true};
                pthread_create(&threads[i], NULL, (void* (*)(void*))bfs_partitions_runner, &bfs_plus_infos[i]);
            }

            for(size_t i = 0; i < threads_to_use; i++) {
                pthread_join(threads[i], NULL);
            }
        }
        SCC_count += unique_colors.size();
        DEB("Finished BFS")

        // remove all vertices that are in some SCC
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
        if (USE_ONB) {
           SCC_count += trimVertices_inplace(inb, onb, vleft, SCC_id, SCC_count);
        } else {
           SCC_count += trimVertices_inplace_single_direction(inb, vleft, SCC_id, SCC_count);
        }
        // clean up vleft after trim
        std::erase_if(vleft, [&](size_t v) { return SCC_id[v] != UNCOMPLETED_SCC_ID; });
    }
    DEB("Finished")

    DEB("Total tries: " << total_tries)
    DEB("Total iterations: " << iter)
    DEB("Total SCCs: " << SCC_count)
    return SCC_id;
}