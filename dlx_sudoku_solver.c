/* this code is a very scrappy (and not well tested) copied and modified version of Knuth's https://www-cs-faculty.stanford.edu/~knuth/programs/dance.w */

#include "stdio.h"
#include "stdlib.h"
#include "stdint.h"
#include "assert.h"

typedef uint8_t boolean;
typedef uint8_t byte;
typedef uint32_t u32;

#define GORNISHT 0
#define MAX_ROWS_IN_COVER 1234

/* all lists are non-empty circular doubly linked lists */

/*
 * each 1 entry in the exact cover matrix gets a node
 * also each column holds a dummy node
 */
struct dlx_node
{
    u32 node_id; /* for debugging */

    u32 row_id; /* for printing the solution */
    struct dlx_node* left; /* left and right pointers won't change as rows won't change fixed */
    struct dlx_node* right;

    struct dlx_node* up;
    struct dlx_node* down;
    struct dlx_col* col;
};

struct dlx_col
{
    struct dlx_node head; /* dummy node as the header of the list of the columns node (available to cover it) */
    u32 len; /* the number of (non-header) nodes */
    u32 col_id; /* for debugging */
    
    struct dlx_col* prev;
    struct dlx_col* next;
};

static u32 make_cntr_mat(boolean* mat, u32* cntr_mat, u32 nrows, u32 ncols)
{
    u32 cntr = 0;
    u32 idx;
    for(idx = 0; idx < nrows * ncols; idx += 1)
    {
        if(mat[idx] != 0)
        {
            assert(mat[idx] == 1);
            cntr += 1;
            cntr_mat[idx] = cntr;
        }
        else { cntr_mat[idx] = 0; }
    }
    return cntr;
}

static void set_dummy_node(struct dlx_node* node)
{
    node->node_id = node->row_id = UINT32_MAX;
    node->left = node->right = node->up = node->down = GORNISHT;
    node->col = GORNISHT;
}

static void set_node_ids(struct dlx_node* node_arr, u32 cntr)
{
    u32 idx;
    for(idx = 1; idx <= cntr; idx += 1) { node_arr[idx].node_id = idx; }
    set_dummy_node(&(node_arr[0]));
}

static void link_rows(u32* cntr_mat, struct dlx_node* node_arr, u32 nrows, u32 ncols, u32* scratch)
{
    u32 num_nodes, node_id, i, j, k;

    for(i = 0; i < nrows; i += 1)
    {
        num_nodes = 0;
        for(j = 0; j < ncols; j += 1)
        {
            node_id = cntr_mat[i * ncols + j];
            if(node_id != 0)
            {
                node_arr[node_id].row_id = i;
                scratch[num_nodes] = node_id;
                num_nodes += 1;
            }
        }
    
        for(k = 0; k < num_nodes; k += 1)
        {
            node_arr[scratch[k]].right = &(node_arr[scratch[(k+1) % num_nodes]]);
            node_arr[scratch[(k+1) % num_nodes]].left = &(node_arr[scratch[k]]);
        }
    }
}

static void link_cols(u32* cntr_mat, struct dlx_node* node_arr, struct dlx_col* col_arr, struct dlx_col* root, u32 nrows, u32 ncols, u32* scratch)
{
    u32 num_nodes, node_id, i, j, k;
    for(j = 0; j < ncols; j += 1)
    {
        col_arr[j].col_id = j;
        set_dummy_node(&(col_arr[j].head));
        col_arr[j].head.col = &(col_arr[j]);

        num_nodes = 0;
        for(i = 0; i < nrows; i += 1)
        {
            node_id = cntr_mat[i * ncols + j];
            if(node_id != 0)
            {
                scratch[num_nodes] = node_id;
                num_nodes += 1;

                node_arr[node_id].col = &(col_arr[j]);
            }
        }
        col_arr[j].len = num_nodes;

        if(num_nodes == 0) { col_arr[j].head.down = col_arr[j].head.up = &(col_arr[j].head); }
        else
        {
            for(k = 0; k + 1 < num_nodes; k += 1)
            {
                node_arr[scratch[k]].down = &(node_arr[scratch[k+1]]);
                node_arr[scratch[k+1]].up = &(node_arr[scratch[k]]);
            }
            col_arr[j].head.down = &(node_arr[scratch[0]]);
            node_arr[scratch[0]].up = &(col_arr[j].head);
            col_arr[j].head.up = &(node_arr[scratch[num_nodes - 1]]);
            node_arr[scratch[num_nodes - 1]].down = &(col_arr[j].head);
        }
    }

    set_dummy_node(&(root->head)); /* we're never supposed to use this node */
    root->len = root->col_id = UINT32_MAX;
    assert(ncols > 0);
    root->next = &(col_arr[0]);
    col_arr[0].prev = root;
    root->prev = &(col_arr[ncols - 1]);
    col_arr[ncols - 1].next = root;

    for(k = 0; k + 1 < ncols; k += 1)
    {
        col_arr[k].next = &(col_arr[k + 1]);
        col_arr[k + 1].prev = &(col_arr[k]);
    }
}

static void cover(struct dlx_col* c)
{

    struct dlx_col* l;
    struct dlx_col* r;
    struct dlx_node* rr;
    struct dlx_node* nn;
    struct dlx_node* uu;
    struct dlx_node* dd;

    assert(c->col_id != UINT32_MAX);

    /* remove column c from uncovered columns */
    l = c->prev;
    r = c->next;
    l->next = r;
    r->prev = l;

    for(rr = c->head.down; rr != &(c->head); rr = rr->down)
    {
        for(nn = rr->right; nn != rr; nn = nn->right)
        {
            /* removing node nn from other columns, as c is already covered */
            uu = nn->up;
            dd = nn->down;
            uu->down = dd;
            dd->up = uu;

            assert(nn->col->len != 0);
            nn->col->len -= 1;
        }
    }
}

static void uncover(struct dlx_col* c)
{
    struct dlx_col* l;
    struct dlx_col* r;
    struct dlx_node* rr;
    struct dlx_node* nn;
    struct dlx_node* uu;
    struct dlx_node* dd;

    assert(c->col_id != UINT32_MAX);

    for(rr = c->head.up; rr != &(c->head); rr = rr->up)
    {
        for(nn = rr->left; nn != rr; nn = nn->left)
        {
            /* insert node nn, it can be used to cover again */
            uu = nn->up;
            dd = nn->down;
            uu->down = nn;
            dd->up = nn;

            nn->col->len += 1;
        }
    }

    /* insert column c into uncovered columns */
    l = c->prev;
    r = c->next;
    l->next = c;
    r->prev = c;
}

static struct dlx_col* find_best_col(struct dlx_col* root)
{
    u32 minlen = UINT32_MAX;
    struct dlx_col* cur_col;
    struct dlx_col* ret = GORNISHT;
    
    assert(root->next != root);

    for(cur_col = root->next; cur_col != root; cur_col = cur_col->next)
    {
        assert(cur_col->col_id != UINT32_MAX);
        if(cur_col->len < minlen)
        {
            ret = cur_col;
            minlen = cur_col->len;
        }
    }
    assert(ret != GORNISHT);
    return ret;
}

static u32 dlx(struct dlx_col* root, boolean stop_at_first_solution, void exact_cover_callback(u32, u32*, byte*), byte* extra)
{
    u32 ret = 0, level = 0, ii; 
    struct dlx_col* best_col;
    struct dlx_node* cur_node;
    struct dlx_node* pp;
    
    struct dlx_node* history[MAX_ROWS_IN_COVER]; /* this array is a stack, its number of elements is level */
    u32 history_row_ids[MAX_ROWS_IN_COVER];

forward:
    best_col = find_best_col(root);
    cover(best_col);
    cur_node = history[level] = best_col->head.down;
    history_row_ids[level] = history[level]->row_id;

advance:
    assert(cur_node->col == best_col);

    if(cur_node == &(best_col->head))
        goto backup;


    for(pp = cur_node->right; pp != cur_node; pp = pp->right)
        cover(pp->col);

    if(root->next == root)
    {
        ret += 1;
        if(exact_cover_callback != GORNISHT)
            exact_cover_callback(level + 1, history_row_ids, extra);
        else
        {
            printf("solution: ");
            for(ii = 0; ii <= level; ii += 1)
                printf("%d ", history_row_ids[ii]);
            printf("\n");
        }
        
        if(stop_at_first_solution) { goto done; }
        goto recover;
    }

    level += 1;

    assert(level < MAX_ROWS_IN_COVER);
    goto forward;

backup:
    uncover(best_col);
    if(level == 0)
        goto done;

    level -= 1;
    cur_node = history[level];
    best_col = cur_node->col;

recover: 
    for(pp = cur_node->left; pp != cur_node; pp = pp->left)
        uncover(pp->col);

    cur_node = history[level] = cur_node->down;
    history_row_ids[level] = history[level]->row_id;

    goto advance;

done:
    return ret;
}


static u32 find_exact_covers(boolean* mat, u32 nrows, u32 ncols, boolean stop_at_first_solution, void exact_cover_callback(u32, u32*, byte*), byte* extra)
{
    assert(nrows > 0);
    assert(ncols > 0);

    u32 ret;
    u32* cntr_mat = malloc(sizeof(u32) * nrows * ncols);
    u32 cntr = make_cntr_mat(mat, cntr_mat, nrows, ncols);
    
    struct dlx_node* node_arr = malloc(sizeof(struct dlx_node) * (cntr + 1)); /* we're being wasetful to avoid off by one errors */
    set_node_ids(node_arr, cntr);
    
    u32 max_dim = nrows > ncols ? nrows : ncols;
    u32* scratch = malloc(sizeof(u32) * max_dim);

    link_rows(cntr_mat, node_arr, nrows, ncols, scratch);

    struct dlx_col* col_arr = malloc(sizeof(struct dlx_col) * ncols);
    struct dlx_col root;
    link_cols(cntr_mat, node_arr, col_arr, &root, nrows, ncols, scratch);
    
    ret = dlx(&root, stop_at_first_solution, exact_cover_callback, extra);

    free(col_arr);
    free(scratch);
    free(node_arr);
    free(cntr_mat); 

    return ret;   
}


static void sudoku_ecm_helper(boolean* ecm, u32 s_col, u32 s_row, u32 s_idx, u32 s_box, u32 s_val)
{
    assert((1 <= s_val) && (s_val <= 9));
    s_val -= 1;

    u32 ecm_row = s_val + 9 * s_col + 9 * 9 * s_row;
    assert(ecm_row < 9 * 9 * 9);

    /* cover cell */
    u32 ecm_col = s_idx;
    assert(ecm_col < 9 * 9);
    ecm[ecm_col + (9 * 9 + 9 * 9 + 9 * 9 + 9 * 9) * ecm_row] = 1;

    /* cover row val */
    ecm_col = (9 * 9) + s_row + 9 * s_val;
    assert(ecm_col - 9 * 9 < 9 * 9);
    ecm[ecm_col + (9 * 9 + 9 * 9 + 9 * 9 + 9 * 9) * ecm_row] = 1;

    /* cover col val */
    ecm_col = (9 * 9) + (9 * 9) + s_col + 9 * s_val;
    assert(ecm_col - 9 * 9 * 2 < 9 * 9);
    ecm[ecm_col + (9 * 9 + 9 * 9 + 9 * 9 + 9 * 9) * ecm_row] = 1;

    /* cover box val */
    ecm_col = (9 * 9) + (9 * 9) + (9 * 9) + s_box + 9 * s_val;
    assert(ecm_col - 9 * 9 * 3 < 9 * 9);
    ecm[ecm_col + (9 * 9 + 9 * 9 + 9 * 9 + 9 * 9) * ecm_row] = 1;
}

static void sudoku_cb(u32 n, u32* arr, byte* extra)
{
    u32* sudoku = (u32*) extra;
    u32 x,i,j,v;
    for(u32 idx = 0; idx < n; idx += 1)
    {
        x = arr[idx];
        v = (x % 9) + 1;
        j = (x / 9) % 9;
        i = x / (9 * 9);
        assert(i < 9);
        assert(j < 9);
        sudoku[j + 9 * i] = v;
    }
    for(i = 0; i < 9; i += 1)
    {
        for(j = 0; j < 9; j += 1)
        {
            printf("%d ", sudoku[j + 9 * i]);
        }
        printf("\n");
    }
    printf("\n");
}

static void sudoku_solver(u32 grid[81])
{
    /*                      val col row     cell  row val col val box val */
    boolean exact_cover_mat[(9 * 9 * 9) * (9 * 9 + 9 * 9 + 9 * 9 + 9 * 9)] = {0};

    u32 s_row, s_col, s_idx, s_box, s_val;

    for(s_row = 0; s_row < 9; s_row += 1)
    {
        for(s_col = 0; s_col < 9; s_col += 1)
        {
            s_idx = s_col + 9 * s_row;
            s_box = (s_col / 3) + 3 * (s_row / 3);
            s_val = grid[s_idx];
            if(s_val != 0) sudoku_ecm_helper(exact_cover_mat, s_col, s_row, s_idx, s_box, s_val);
            else
            {
                for(s_val = 1; s_val <= 9; s_val += 1)
                {
                    sudoku_ecm_helper(exact_cover_mat, s_col, s_row, s_idx, s_box, s_val);
                }
            }
        }
    }

    find_exact_covers(exact_cover_mat, 9 * 9 * 9, 9 * 9 + 9 * 9 + 9 * 9 + 9 * 9, 0, sudoku_cb, (byte*) grid);
}

int main()
{
    u32 grid[81] = {
        2,0,0,0,0,5,0,0,9,
        0,7,0,4,0,0,5,0,6,
        0,0,6,1,0,0,0,0,0,
        5,0,0,0,8,0,0,0,2,
        0,0,3,0,0,0,1,0,0,
        1,0,0,0,9,0,0,0,7,
        0,0,0,0,0,2,9,0,0,
        4,0,2,0,0,3,0,7,0,
        8,0,0,9,0,0,0,0,1};

    sudoku_solver(grid);

    return 0;
}
