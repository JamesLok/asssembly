#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <assert.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"
#include "cmd_line.h"
#include "graph_info.h"
#include "db_differentiation.h"
#include "db_complex_genotyping.h"
#include "model_selection.h"
#include "experiment.h"
#include "genome_complexity.h"
#include "maths.h"
#include "seq_error_rate_estimation.h"

// third party libraries
#include <seq_file.h>
#include <string_buffer.h>

#include "binary_kmer.h"

#define QN 400
#define NN 16000
#define NN_PATH 1000

#define mode 0

//typedef struct Node_with_Orient
//{
//  dBNode* node;
//  Orientation orient;
//  Node_with_Orient *previous;
//}Node_with_Orientation;

typedef struct
{
  dBNode* node;
  Orientation orient;
  dBNode *pre_node;
  Orientation pre_orient;
}Node_with_Orientation;


typedef struct Queue  
{  
    Node_with_Orientation * pBase;    
    int front;       
    int rear;         
}QUEUE;


boolean isFull(QUEUE * pQueue)  
{  
    if((pQueue->rear+1) % QN == pQueue->front)     //队列满  
        return true;  
    else  
        return false;  
}  
/* 
 *判断队列是否为空函数的实现 
 */  
boolean isEmpty(QUEUE * pQueue)  
{  
    if(pQueue->front == pQueue->rear)  
        return true;  
    else  
        return false;  
}

void initQueue(QUEUE * pQueue)  
{  
    //分配内存  
    pQueue->pBase = (Node_with_Orientation *)malloc(sizeof(Node_with_Orientation) * QN);          //分配10个int型所占的空间  
    pQueue->front = 0;       //初始化时，front和rear值均为0  
    pQueue->rear = 0;  
    return;  
}  
/* 
 *入队函数的实现 
 */  
boolean enQueue(QUEUE * pQueue,  Node_with_Orientation *pvalue)  
{  
    if(isFull(pQueue))  
    {  
        printf("队列已满,不能再插入元素了!\n");  
        return false;  
    }  
    else  
    {  
        //向队列中添加新元素  
        pQueue->pBase[pQueue->rear].node = pvalue->node;
        pQueue->pBase[pQueue->rear].orient = pvalue->orient;
        pQueue->pBase[pQueue->rear].pre_node = pvalue->pre_node; 
        pQueue->pBase[pQueue->rear].pre_orient = pvalue->pre_orient; 
        //将rear赋予新的合适的值  
        pQueue->rear = (pQueue->rear+1) % QN;  
        return true;  
    }  
}  
/* 
 *出队函数的实现 
 */  
boolean outQueue(QUEUE * pQueue, Node_with_Orientation *pvalue)  
{  
    //如果队列为空,则返回false  
    if(isEmpty(pQueue))  
    {  
        printf("队列为空，出队失败!\n");  
        return false;  
    }  
    else  
    {  
        pvalue->node = pQueue->pBase[pQueue->front].node;
        pvalue->orient = pQueue->pBase[pQueue->front].orient;
        pvalue->pre_node = pQueue->pBase[pQueue->front].pre_node;
        pvalue->pre_orient = pQueue->pBase[pQueue->front].pre_orient;
        
        pQueue->front = (pQueue->front+1) % QN;      //移到下一位置  
        return true;  
    }  
}


boolean cheack_whether_edge_already_in(Node_with_Orientation *curr_from, Node_with_Orientation *curr_to, Node_with_Orientation *curr_within,
                                       Node_with_Orientation *from_array, Node_with_Orientation *to_array, Node_with_Orientation *within_array, int jNN)
{
    int k;
    for (k=0; k<jNN; k++)
    {
        if (curr_from->node == from_array[k].node && curr_from->orient == from_array[k].orient &&
            curr_to->node == to_array[k].node && curr_to->orient == to_array[k].orient &&
            curr_within->node == within_array[k].node && curr_within->orient == within_array[k].orient ) //supernode already done
        {
            printf("node bf supernode already in.\n");
            return true;
        }
    }
    return false;
}

Nucleotide get_one_edge_all_colours(Element *tmp_node, Orientation orientation)
{
    Edges edge;
    edge = element_get_colour_union_of_all_colours(tmp_node);
    if (orientation == reverse)
    {
       edge >>= 4;
    }
    edge = edge & 0x0F;
    
    if (edge==0) return 4;
    
    if (edge != 8 && edge != 4 && edge != 2 && edge !=1)
        die("more edges exists!\n");
    
    int n;
    for(n=0;n<4;n++)
    {
       if ((edge & 1) == 1)
       {
            return(n);
       }
       edge >>= 1; 
    }
    
    return 4;
}

int count_num_of_edges_all_colours(Element *tmp_node, Orientation orientation)
{
    Edges edge;
    edge = element_get_colour_union_of_all_colours(tmp_node);
    if (orientation == reverse)
    {
       edge >>= 4;
    }
    edge = edge & 0x0F;
    
    int i, num=0;
    char tmp = 1;
    for (i=0; i<4; i++)
    {
        if (edge & tmp)
            num++;
        tmp = tmp <<1;
    }
    return num;
}


Edges db_node_get_edge_in_specific_person_or_population(dBNode * element, Orientation orientation, int colour)
{
  //get the edge char for this specific person or pop:
  char edge = get_edge_copy(*element, colour);

  if (orientation == reverse)
  {
    edge >>= 4;
  }
  
  return (edge & 0x0F);
}

void print_node_to_fasta(char *out_file_name, Node_with_Orientation *node_w_orient, int node_array_len, int kmer_size)
{
    FILE *fout = fopen(out_file_name, "w");
    if(fout == NULL)
    {
      die("cannot open file:%s\n", out_file_name);
    }
    
    int i;
    printf("output node");
    char str[kmer_size+1];
    BinaryKmer tmp_kmer;
    
    for(i=0; i<node_array_len; i++)
    {
        fprintf(fout,">node%d\n",i);
        
        if(node_w_orient[i].orient == forward)
        {
            
            binary_kmer_to_seq((BinaryKmer*)element_get_kmer(node_w_orient[i].node), kmer_size, str);
            str[kmer_size] = '\0';
            fprintf(fout, "%s\n", str);
        }
        else if(node_w_orient[i].orient == reverse)
        {
            binary_kmer_reverse_complement(element_get_kmer(node_w_orient[i].node), kmer_size, &tmp_kmer);
            binary_kmer_to_seq(tmp_kmer, kmer_size, str);
            str[kmer_size] = '\0';
            fprintf(fout,"%s\n", str);
        }
    }
    fclose(fout); 
    
}

void _print_kmer_with_orientation(dBNode* bkmer, short kmer_size, Orientation orientation)
{
    BinaryKmer tmp_kmer;
    if (orientation == forward)
        _print_kmer(element_get_kmer(bkmer), kmer_size);
    else if(orientation == reverse)
    {
        printf("reverse!!!\n");
        _print_kmer(binary_kmer_reverse_complement(element_get_kmer(bkmer), kmer_size, &tmp_kmer ), kmer_size);
        
    }
    printf("\n");
}

                            
Edges db_node_which_edge_exist(dBNode * element,
                           Orientation orientation, int colour)
{
  //get the edge char for this specific person or pop:
  char edge = get_edge_copy(*element, colour);

  if (orientation == reverse)
  {
    edge >>= 4;
  }
  
  return (edge & 0x0F);
}


double get_avg_covg_for_specific_color_from_to_node(Element *from_node, Orientation from_orient, Element *to_node, Orientation to_orient, Element *within_node, Orientation within_orient,
                                                    int max_expected_sup_len, dBGraph* db_graph, int index)
{
    double sum_coverage = 0.0;
    
    if (to_node == within_node && to_orient==within_orient)
        return(0.0);
    

    
    dBNode* curr_node = within_node;
    Orientation curr_orientation = within_orient;
    boolean found = false;
    
    void check_edge(Nucleotide base)
    {
        dBNode * tmp_next_node;
        Orientation tmp_next_orientation;
        Nucleotide rev_base;
        
        if (db_node_edge_exist(from_node, base, from_orient, index)){
                                    
            tmp_next_node = db_graph_get_next_node_for_specific_person_or_pop(from_node,from_orient,&tmp_next_orientation,
                                                                                                base,&rev_base,db_graph, index);
            printf("color = %d, next_node:  ", index);
            _print_kmer_with_orientation(tmp_next_node,19, tmp_next_orientation);
            if (tmp_next_node == curr_node && tmp_next_orientation == curr_orientation)
                found = true;
        }                 
    }
    nucleotide_iterator(&check_edge);
    
    if (found == false)
    {
        return(0.0);
    }
    
    
    int super_lenth = 0;
    dBNode* next_node;
    Nucleotide nuc, reverse_nuc;   
    Orientation next_orientation;
    do
    {
        
         if(!db_node_has_precisely_one_edge(curr_node, curr_orientation, &nuc, index))
            die("cal_covg: node_has_more_edges\n");
         sum_coverage = sum_coverage + db_node_get_coverage(curr_node, index);
         super_lenth++;
         
         next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orientation, &next_orientation, nuc, &reverse_nuc,db_graph, index);
         
         curr_node = next_node;
         curr_orientation = next_orientation;
         
//        _print_kmer_with_orientation(curr_node,31, curr_orientation);
//        _print_kmer_with_orientation(to_node,31, to_orient);
          
         
        
    }while((next_node != to_node) || (next_orientation != to_orient));
    
    if (super_lenth>max_expected_sup_len)
    {
        die("cal_covg: len > max_expected_sup_len\n");
    }
    printf("len = %d\n", super_lenth);
    return(sum_coverage/super_lenth);
    
}


boolean get_supernode_avg_covg_for_each_color(Element *curr_node, Orientation curr_orient, int max_expected_sup_len,
                                               dBGraph* db_graph, double *mat, int offset)
{
    dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)* (max_expected_sup_len+1));
    Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*(max_expected_sup_len+1)); 
    Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_expected_sup_len+1));
    char*        supernode_string  = (char*) malloc(sizeof(char)*(max_expected_sup_len+2)); //2nd +1 for \0
    double avg_cov = 0.0;
    Covg min_cov;
    Covg max_cov;
    boolean is_cycle;
    
    
    int i;
    
    for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
        db_graph_supernode_for_specific_person_or_pop(curr_node,max_expected_sup_len,
                                                    &db_node_action_do_nothing,
                                                    path_nodes, path_orientations, path_labels, supernode_string,
                                                    &avg_cov,&min_cov, &max_cov, &is_cycle,
                                                    db_graph, i);
        printf("%f\n",avg_cov);
        *(mat + offset + i) = avg_cov;
       
    }
    
    
    free(path_nodes);
    free(path_orientations);
    free(path_labels);
    free(supernode_string);
    
    return true;
}



boolean get_one_supernode_and_prepare_for_next(Element **pcurr_node, Orientation *pcurr_orient, int max_expected_sup_len,
                                               dBGraph* db_graph, Node_with_Orientation *pfrom, Node_with_Orientation *pto, Node_with_Orientation *pwithin, int *psuper_len,
                                               Node_with_Orientation *node_from, Node_with_Orientation *node_to, Node_with_Orientation *node_within_super, int *supernode_len,  int *current_jNN)
{
    dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)* (max_expected_sup_len+1));
    Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*(max_expected_sup_len+1)); 
    Nucleotide*  path_labels       = (Nucleotide*) malloc(sizeof(Nucleotide)*(max_expected_sup_len+1));
    char*        supernode_string  = (char*) malloc(sizeof(char)*(max_expected_sup_len+2)); //2nd +1 for \0
    double avg_cov = 0.0;
    Covg min_cov;
    Covg max_cov;
    boolean is_cycle;
    int i;

    
    Orientation curr_orient = *pcurr_orient;
    Element *curr_node = *pcurr_node;
    
//                    _print_kmer_with_orientation(curr_node,kmer_size, curr_orient);
    
    
    if (db_node_check_status(curr_node, visited)==true)
    {
        free(path_nodes);
        free(path_orientations);
        free(path_labels);
        free(supernode_string);
        
//        *psuper_len = -1;
        
        return false;
    }
    
    int length_sup =  db_graph_supernode_in_subgraph_defined_by_func_of_colours(curr_node,max_expected_sup_len,
                                                            &db_node_action_set_status_visited,
                                                            path_nodes, path_orientations, path_labels, supernode_string,
                                                            &avg_cov,&min_cov, &max_cov, &is_cycle,
                                                            db_graph, 
                                                            &element_get_colour_union_of_all_colours, &element_get_covg_union_of_all_covgs);


    if (curr_orient == reverse)
    {
        dBNode* tmp;
        Orientation tmp_orient1, tmp_orient2;
        
        //tmp = path_nodes[0];
        //path_nodes[0] = path_nodes[length_sup];
        //path_nodes[length_sup] = tmp;
        //
        //Orientation tmp_orient_1st, tmp_orient_last;
        //tmp_orient_1st = path_orientations[0];
        //tmp_orient_last = path_orientations[length_sup];
        //if (tmp_orient_1st == forward)
        //    path_orientations[length_sup] = reverse;
        //else if (tmp_orient_1st == reverse)
        //    path_orientations[length_sup] = forward;
        //
        //if (tmp_orient_last == forward)
        //    path_orientations[0] = reverse;
        //else if (tmp_orient_last == reverse)
        //    path_orientations[0] = forward;
        

        for(i=0; i< (int)(length_sup/2 +1);i++)
        {
            tmp = path_nodes[i];
            path_nodes[i] = path_nodes[length_sup-i];
            path_nodes[length_sup-i] = tmp;
        }
        
        for(i=0; i< (int)(length_sup/2 +1);i++)
        {
            tmp_orient1 = path_orientations[i];
            tmp_orient2 = path_orientations[length_sup-i];
            if (tmp_orient1 == forward)
                path_orientations[length_sup-i] = reverse;
            else if (tmp_orient1 == reverse)
                path_orientations[length_sup-i] = forward;
            
            if (tmp_orient2 == forward)
                path_orientations[i] = reverse;
            else if (tmp_orient2 == reverse)
                path_orientations[i] = forward;
        }
    }
    
    
    //need to check whehter the color profile along the path are the same
    
    boolean color_edge_match = true, path_segmentation = false;
    
    if (length_sup>1)
    {
    
        int colour;
        Edges edge_colour[NUMBER_OF_COLOURS];
            
        for (colour=0; colour<NUMBER_OF_COLOURS;colour++)
            edge_colour[colour] = db_node_which_edge_exist(path_nodes[1],path_orientations[1],colour)>0 ? 1 : 0;
        
        
        int start = 0;
        
        for (i=1; i<length_sup;i++)
        {
            for (colour=0; colour<NUMBER_OF_COLOURS;colour++)
            {
    //            printf("colour (%d) = %d\n",colour,db_node_which_edge_exist(path_nodes[i],path_orientations[i],colour));
     //           printf("colour (%d) = %d\n",colour,db_node_which_edge_exist(path_nodes[i],path_orientations[i],colour)>0 ? 1 : 0);
                if (edge_colour[colour] != (db_node_which_edge_exist(path_nodes[i],path_orientations[i],colour)>0 ? 1 : 0))
                {
                    color_edge_match = false;
                    path_segmentation = true;
                    break;
                    
                }
                
            }
    //        printf("%d\n", *current_jNN);
    //        printf("colour (%d) = %d\n",0,db_node_which_edge_exist(path_nodes[i],path_orientations[i],0));
    //       printf("colour (%d) = %d\n",1,db_node_which_edge_exist(path_nodes[i],path_orientations[i],1));
            if (!color_edge_match)
            {
                printf("in it\n");
                printf("%d (%d, %d | %d)\n", (*current_jNN), start,i, length_sup);
                node_from[*current_jNN].node = path_nodes[start];
                node_from[*current_jNN].orient = path_orientations[start];
                
                node_to[*current_jNN].node = path_nodes[i];
                node_to[*current_jNN].orient = path_orientations[i];
                
    //            dBNode * tmp_next_node;
    //            Orientation tmp_next_orientation;
    //            Nucleotide rev_base;
              
    //            Nucleotide base = get_one_edge_all_colours(path_nodes[start], path_orientations[start]);
    
    //            tmp_next_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(path_nodes[start], path_orientations[start],&tmp_next_orientation,
    //                                                                                    base,&rev_base,db_graph, &element_get_colour_union_of_all_colours);
                
                node_within_super[*current_jNN].node = path_nodes[start+1];
                node_within_super[*current_jNN].orient = path_orientations[start+1];
                *(supernode_len + (*current_jNN)) = i - start - 1;
                
                db_node_set_status_to_none(node_from[*current_jNN].node);
                db_node_set_status_to_none(node_to[*current_jNN].node);
                
                boolean pre_from_to_node_found = cheack_whether_edge_already_in(&node_from[*current_jNN], &node_to[*current_jNN], &node_within_super[*current_jNN], node_from, node_to, node_within_super,*current_jNN);
                                
                                
                if (!pre_from_to_node_found)            
                    (*current_jNN) = (*current_jNN) + 1;
                
                
                start = i;
                
                for (colour=0; colour<NUMBER_OF_COLOURS;colour++)
                    edge_colour[colour] = db_node_which_edge_exist(path_nodes[i],path_orientations[i],colour)>0 ? 1 : 0;
            
                color_edge_match = true;
                
                
            }
     
            
        }
        
        if (color_edge_match & path_segmentation)
        {
            printf("%d (%d, %d | %d)\n", (*current_jNN), start,i, length_sup);
            
            node_from[*current_jNN].node = path_nodes[start];
            node_from[*current_jNN].orient = path_orientations[start];
            
            node_to[*current_jNN].node = path_nodes[length_sup];
            node_to[*current_jNN].orient = path_orientations[length_sup];
            
            //dBNode * tmp_next_node;
            //Orientation tmp_next_orientation;
            //Nucleotide rev_base;
            //
            //Nucleotide base = get_one_edge_all_colours(path_nodes[start], path_orientations[start]);
            //
            //tmp_next_node = db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(path_nodes[start], path_orientations[start],&tmp_next_orientation,
            //                                                                        base,&rev_base,db_graph, &element_get_colour_union_of_all_colours);
            
            node_within_super[*current_jNN].node = path_nodes[start+1];
            node_within_super[*current_jNN].orient = path_orientations[start+1];
            *(supernode_len + (*current_jNN)) = length_sup - start - 1;
            
            
            db_node_set_status_to_none(node_from[*current_jNN].node);
            db_node_set_status_to_none(node_to[*current_jNN].node);
            
            boolean from_to_node_found = cheack_whether_edge_already_in(&node_from[*current_jNN], &node_to[*current_jNN], &node_within_super[*current_jNN], node_from, node_to, node_within_super,*current_jNN);
                            
                            
            if (!from_to_node_found)            
                (*current_jNN) = (*current_jNN) + 1;      
        }
        
    }
    
    printf("%d\n", length_sup);
    printf("%s\n", supernode_string);
    //printf("\n-------------------------------------\n");
    //
    //_print_kmer_with_orientation(path_nodes[0],kmer_size, path_orientations[0]);
    //
    //
    //printf("\n-------------------------------------\n");
    //_print_kmer_with_orientation(path_nodes[length_sup-1],kmer_size, path_orientations[length_sup-1]);
    //printf("\n-------------------------------------\n");
    //
    _print_kmer_with_orientation(path_nodes[length_sup],19, path_orientations[length_sup]);
    printf("\n---------------%f----------------------\n",avg_cov);
    
    
    *pcurr_node = path_nodes[length_sup];
    *pcurr_orient = path_orientations[length_sup];

    pfrom->node = path_nodes[0];
    pfrom->orient = path_orientations[0];
    
    pto->node = path_nodes[length_sup];
    pto->orient = path_orientations[length_sup];
    
    pwithin->node = path_nodes[(length_sup>0)?1:0];
    pwithin->orient = path_orientations[(length_sup>0)?1:0];
    
    *psuper_len = length_sup-1;
    
    /////////////////////////////////////
    // test
    if (length_sup == 1)
    {
        printf("length_sup == 1!!!!\n");
        _print_kmer_with_orientation(path_nodes[0],19, path_orientations[0]);
    }
    
    ////////////////////////////////
    
    
    free(path_nodes);
    free(path_orientations);
    free(path_labels);
    free(supernode_string);
    
    if (path_segmentation)
    {
        printf("path_segmentation!\n");
        return false;
    }   
    
    return true;
}


void reduce_bubble_to_one_path(dBGraph* db_graph_small, Orientation tip_orient)
{
    int i;
    int count_forward, count_reverse;
    int next_node_in, next_node_out;
}

void remove_shorter_tips(dBGraph* db_graph_small, Orientation tip_orient)
{
    int i;
    int count_forward, count_reverse;
    int next_node_in, next_node_out;
    
    dBNode *curr_node, *next_node;
    Orientation curr_orient, next_orientation;
    Nucleotide nuc, reverse_nuc;
    Edges edge;
    int n;
    
    
    for(i=0;i < db_graph_small->number_buckets * db_graph_small->bucket_size;i++)
    {
        if (!db_node_check_for_flag_ALL_OFF(&db_graph_small->table[i]) && db_node_check_status(&db_graph_small->table[i], none))
        {
            
            count_forward = count_num_of_edges_all_colours(&db_graph_small->table[i],forward);
            count_reverse = count_num_of_edges_all_colours(&db_graph_small->table[i],reverse);    
            
            if ((count_forward == 0 && count_reverse == 1 && tip_orient == reverse) || (count_forward == 1 && count_reverse == 0 && tip_orient == forward))
            {
                int total_len = 0;
                int total_covg = 0;
                curr_node = &db_graph_small->table[i];
                curr_orient = tip_orient;
                printf("staring a new one = %d!\n", db_node_get_coverage(curr_node, 0));
                do
                {
                    db_node_set_status(curr_node, visited);
                    total_covg += db_node_get_coverage(curr_node, 0);
                    
                    edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                    
                    for(n=0;n<4;n++)
                    {
                        if ((edge & 1) == 1)
                        {
                            nuc = n;

                            next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                            if (next_node == NULL)
                                die("somthing wrong in small db_graph(1)\n");
                           
                        }
                        
                        edge >>= 1;    
                    }
                    
                    //printf("%d\n", total_len);
                    total_len++;
                    //printf("node %p status = %d\n", curr_node, curr_node->status);
                    //_print_kmer_with_orientation(curr_node, kmer_size, curr_orient);
                    curr_node = next_node;
                    curr_orient = next_orientation;
                    
                    
                }while(count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 1 && count_num_of_edges_all_colours(curr_node,curr_orient) == 1);
                
                printf("finish one branch! total_len = %d\n", total_len);
                printf("node coverage = %d\n", db_node_get_coverage(curr_node, 0));
                printf("node_in = %d, node_out = %d\n", count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)), count_num_of_edges_all_colours(curr_node,curr_orient));
                //_print_kmer_with_orientation(curr_node, kmer_size, curr_orient);
                if (count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 2)
                {
                    curr_orient = opposite_orientation(curr_orient);
                    
                    dBNode *junction_node  = curr_node;
                    Orientation junction_orient = curr_orient;
                    
                    int junction_nuc[2];
                    
                    edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                    printf("edge = %d\n", edge);
                    for(n=0;n<4;n++)
                    {
                        if ((edge & 1) == 1)
                        {
                            nuc = n;

                            next_node = db_graph_get_next_node_for_specific_person_or_pop(junction_node, junction_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                            //printf("next node %p status = %d\n", next_node, next_node->status);
                            //_print_kmer_with_orientation(next_node, kmer_size, next_orientation);
                            if (next_node == NULL)
                                {printf("total_len = %d\n", total_len); die("somthing wrong in small db_graph(2)\n");}
                                
                            if (!db_node_check_status(next_node, visited))
                            {
                                curr_node = next_node;
                                curr_orient = next_orientation;
                                junction_nuc[1] = nuc;
                            //    printf("not visited, nuc = %d", nuc);
                                
                            }
                            else if (db_node_check_status(next_node, visited))
                            {
                                junction_nuc[0] = nuc;
                            //    printf("junc_nuc = %d\n",nuc);
                            }
                           
                        }
                        
                        edge >>= 1;    
                    }
                    
                    printf("finish finding the another direction, node = %d\n", db_node_get_coverage(curr_node, 0));
                    // go to another branch
                    do
                    {
                        total_covg -= db_node_get_coverage(curr_node, 0);
                        db_node_set_status(curr_node, visited);
                        
                        edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                        
                        for(n=0;n<4;n++)
                        {
                            if ((edge & 1) == 1)
                            {
                                nuc = n;
    
                                next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                                if (next_node == NULL)
                                    die("somthing wrong in small db_graph(3)\n");
                               
                            }
                            
                            edge >>= 1;    
                        }
                        
                        //printf("%d\n", total_len);
                        total_len--;
                        curr_node = next_node;
                        curr_orient = next_orientation;
                        

                    }while(count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 1 && count_num_of_edges_all_colours(curr_node,curr_orient) == 1);

                    
                    printf("finish travelling the other path, node = %d, total_len = %d\n", db_node_get_coverage(curr_node, 0), total_len);
                    if (total_len < 0)
                    {
                        nuc = junction_nuc[0];
                    }
                    else if ((total_len > 0) && (count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 1 && count_num_of_edges_all_colours(curr_node,curr_orient) == 0))
                    {
                        nuc = junction_nuc[1];
                    }
                    else if ((total_len == 0) && (count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 1 && count_num_of_edges_all_colours(curr_node,curr_orient) == 0))
                    {
                        if (total_covg >= 0)
                        {
                            nuc = junction_nuc[1];
                        }
                        else if (total_covg < 0)
                        {
                            nuc = junction_nuc[0];
                        }
                    }
                    else
                        printf("wrong somewhere\n");
                        
                    //printf("junc_nuc = %d\n",nuc);
                    next_node = db_graph_get_next_node_for_specific_person_or_pop(junction_node, junction_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                    //printf("next_node = %p, total_len = %d\n", next_node, total_len);
                    
                    db_node_reset_edge(junction_node,junction_orient, nuc, 0);
                    db_node_reset_edge(next_node,opposite_orientation(next_orientation), reverse_nuc, 0);
                    
                    curr_node = next_node;
                    curr_orient = next_orientation;
                    
                    printf("finish reseting the junction node\n");
                    
                    do
                    {
                        db_node_set_status(curr_node, pruned);
                        
                        edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                        
                        for(n=0;n<4;n++)
                        {
                            if ((edge & 1) == 1)
                            {
                                nuc = n;
    
                                next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                                if (next_node == NULL)
                                    die("somthing wrong in small db_graph(4)\n");
                               
                            }
                            
                            edge >>= 1;    
                        }
                        
                        next_node_in = count_num_of_edges_all_colours(next_node,opposite_orientation(next_orientation));
                        next_node_out = count_num_of_edges_all_colours(next_node,next_orientation);
                        
                        db_node_reset_edge(curr_node,curr_orient, nuc, 0);
                        db_node_reset_edge(next_node,opposite_orientation(next_orientation), reverse_nuc, 0);
                        
                        curr_node = next_node;
                        curr_orient = next_orientation;

                    }while(next_node_in == 1 && next_node_out == 1);
                    
                    printf("finish pruned a tip!\n");
                    
                    if (next_node_out == 0)
                        db_node_set_status(curr_node, pruned);

                    
                }
                printf("%d\n", total_len);
                db_graph_health_check(true, db_graph_small);
                hash_table_traverse(&db_node_action_set_status_of_unpruned_to_none, db_graph_small);
//                if (total_len == -4366)
//                {
//                    if (cmd_line->dump_binary==true)
//                        db_graph_dump_single_colour_binary_of_colour0(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,
//							db_graph_small, db_graph_info, BINVERSION);    
//                }
              
            }
        }
    }
        
}


Covg walk_through_graph(dBGraph* db_graph, dBNode * start_node, Orientation start_orient, int color, boolean *cycle_detected, int *contig_length)
{
    dBNode * curr_node = start_node;
    Orientation curr_orient = start_orient;
    
    Edges edge;
    Covg max_covg = 0, temp_max_covg;
    int n;
    

    Covg sum_coverage = 0;
    
    
    dBNode* next_node;
    dBNode *temp_max_covg_node;
    Orientation next_orientation, temp_next_orient;
    Nucleotide nuc, reverse_nuc;
   
    *cycle_detected = false;
    *contig_length = 0;
    dBNode *next_node_in_cycle = NULL, *cycle_out_node = NULL;
    Orientation next_orient_in_cycle, cycle_out_orient;
    
    
    dBNode * current_node_in_small_graph  = NULL;
    BinaryKmer tmp_kmer, next_kmer;
    
    Covg curr_node_covg;
        
    while ((edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, color)) && db_node_check_status(curr_node,none) )
    {
        curr_node_covg = db_node_get_coverage(curr_node, color);
        
        if (curr_node_covg < 10)
            break;
        
        (*contig_length) = (*contig_length) + 1;
        sum_coverage += curr_node_covg;


        
        max_covg = 0;
        temp_max_covg_node = NULL;
        for(n=0;n<4;n++)
        {
            if ((edge & 1) == 1)
            {
                nuc = n;
              
                next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc, &reverse_nuc, db_graph, color);
                if (next_node == NULL)
                    continue;
                temp_max_covg = db_node_get_coverage(next_node, color);
                if (temp_max_covg > max_covg)
                {
                    temp_max_covg_node = next_node;
                    max_covg = temp_max_covg;
                    temp_next_orient = next_orientation;
                }
            }
            
            edge >>= 1;    
        }
        
        
        db_node_set_status(curr_node,visited);
        
        if (temp_max_covg_node == NULL)
            break;
        //if a cycle is detected
        if (db_node_check_status(temp_max_covg_node, visited))
        {
            printf("in cycle\n");
            
            *cycle_detected = true;
            
            curr_node = temp_max_covg_node;
            curr_orient = temp_next_orient;
            
            max_covg = 0;
            printf("curr_node = %p\n", curr_node);
            _print_kmer_with_orientation(curr_node, db_graph->kmer_size, curr_orient);
            while (edge = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, color) && db_node_check_status(curr_node,visited))
            {
                printf("1\n");
                for(n=0;n<4;n++)
                {
                    if ((edge & 1) == 1)
                    {
                        nuc = n;
                      
                        next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc, &reverse_nuc, db_graph, color);
                        printf("next_node = %p\n", next_node);
                        _print_kmer_with_orientation(next_node, db_graph->kmer_size, next_orientation);
                        if (db_node_check_status(next_node, visited) || db_node_check_status(next_node, special_visited))
                        {
                            next_node_in_cycle = next_node;
                            next_orient_in_cycle = next_orientation;
                            
                            continue;
                        }
                        
                        temp_max_covg = db_node_get_coverage(next_node, color);
                        if (temp_max_covg > max_covg)
                        {
                            temp_max_covg_node = next_node;
                            max_covg = temp_max_covg;
                            temp_next_orient = next_orientation;
                            cycle_out_node = curr_node;
                            cycle_out_orient = curr_orient;
                        }
                        
                        
                        
                    }
                    
                    edge >>= 1;    
                }
                printf("2\n");
                printf("next_node_in_cycle = %p\n", next_node_in_cycle);
                assert(db_node_check_status(next_node_in_cycle, visited) || db_node_check_status(next_node_in_cycle, special_visited));
                printf("3\n");
                db_node_set_status(curr_node,special_visited);
                printf("4\n");
                if (db_node_check_status(next_node_in_cycle, special_visited))
                {
                    curr_node = cycle_out_node;
                    curr_orient = cycle_out_orient;
                    break;
                }
                
                curr_node = next_node_in_cycle;
                curr_orient = next_orient_in_cycle;
                
            }
            
        }
                

        
        //current_node_in_small_graph = hash_table_insert(element_get_key(element_get_kmer(curr_node),db_graph_small->kmer_size, &tmp_kmer),db_graph_small);
        //
        //
        //binary_kmer_assignment_operator(next_kmer, temp_max_covg_node->kmer);
        //
        //
        //if (temp_next_orient == reverse){
        //   binary_kmer_assignment_operator(next_kmer, *(binary_kmer_reverse_complement(&next_kmer,db_graph_small->kmer_size, &tmp_kmer)));
        //}
        //
        //db_node_add_labeled_edge(current_node_in_small_graph, curr_orient, binary_kmer_get_last_nucleotide(&next_kmer), color);
        //
        //
        //if (db_node_check_status(curr_node, visited))
        //    db_node_set_status(current_node_in_small_graph, none);
        //else if (db_node_check_status(curr_node, special_visited))
        //    db_node_set_status(current_node_in_small_graph, special_none);
        
        curr_node = temp_max_covg_node;
        curr_orient = temp_next_orient;
        
        //printf("%c\t%d\t%d\t%d\n", binary_nucleotide_to_char(binary_kmer_get_last_nucleotide(&next_kmer)), max_covg, *contig_length, sum_coverage);
//        _print_kmer_with_orientation(curr_node, db_graph->kmer_size, curr_orient);
    }
    
    
    return (sum_coverage);
    
    

}


int main(int argc, char **argv)
{
    int i;
    
    CmdLine* cmd_line = cmd_line_alloc();
    if (cmd_line==NULL)
      {
        die("Out of memory!! Cannot even malloc the space to store your command-line arguments. Check who else is using your server, you seem to have severe problems\n");
      }
    
    parse_cmdline(cmd_line, argc,argv,sizeof(Element));
  
    int hash_key_bits, bucket_size;
    dBGraph * db_graph = NULL;
    short kmer_size;

    //set hash table variables:
    kmer_size        = cmd_line->kmer_size;
    hash_key_bits    = cmd_line->number_of_buckets_bits; //number of buckets: 2^hash_key_bits
    bucket_size      = cmd_line->bucket_size;
  
    int number_of_bitfields = ((kmer_size * 2) / (sizeof(bitfield_of_64bits)*8))+1;
    int max_kmer_size = (NUMBER_OF_BITFIELDS_IN_BINARY_KMER * sizeof(bitfield_of_64bits) * 4) -1;
    int min_kmer_size = ((NUMBER_OF_BITFIELDS_IN_BINARY_KMER-1) * sizeof(bitfield_of_64bits) * 4) + 1;
  
    if (number_of_bitfields != NUMBER_OF_BITFIELDS_IN_BINARY_KMER)
    {
      die("K-mer %i  is not in current range of kmers [%i - %i] required for this executable!\n",
          kmer_size, min_kmer_size, max_kmer_size);
    }
    
    printf("Maximum k-mer size (compile-time setting): %i\n", max_kmer_size);
  
    if(cmd_line->kmer_size > max_kmer_size)
    {
      die("k-mer size is too big [%i]!",cmd_line->kmer_size);
    }
  
    printf("Actual K-mer size: %d\n", cmd_line->kmer_size);
    
    int max_expected_sup_len = cmd_line->max_var_len;
    
    //Create the de Bruijn graph/hash table
    int max_retries=15;
    db_graph = hash_table_new(hash_key_bits,bucket_size, max_retries, kmer_size);
    if (db_graph==NULL)
    {
        die("Giving up - unable to allocate memory for the hash table\n");
    }
    printf("Hash table created, number of buckets: %d\n",1 << hash_key_bits);
  
  
  
    GraphInfo* db_graph_info=graph_info_alloc_and_init();//will exit it fails to alloc.

    int first_colour_data_starts_going_into=0;

    if (cmd_line->input_multicol_bin==true)
    {
        long long  bp_loaded = load_multicolour_binary_from_filename_into_graph(cmd_line->multicolour_bin,db_graph, 
                                                                                db_graph_info, &first_colour_data_starts_going_into);
        printf("Loaded the multicolour binary %s, and got %qd kmers\n", cmd_line->multicolour_bin, bp_loaded/db_graph->kmer_size);
    }
    
    
    
    // find the k-mer with the largest coverage
    Covg temp_max, covg_max = 0;
    int index_max_covg = 0;
    
    int num_small_contig = 0;
    dBNode * curr_node ;
    Orientation curr_orient;
    
    while(true)
    {
        covg_max = 0;
        for(i=0;i < db_graph->number_buckets * db_graph->bucket_size;i++)
        {
          if (!db_node_check_for_flag_ALL_OFF(&db_graph->table[i]) && db_node_check_status(&db_graph->table[i],none))
          {
            
            temp_max = db_node_get_coverage(&db_graph->table[i], 0);
            if (covg_max < temp_max)
            {
              covg_max = temp_max;
              index_max_covg = i;
            }
          }  
        }
        printf("finish search\n");
        curr_node = &db_graph->table[index_max_covg];
        
        
        
        int contig_length_forward, contig_length_reverse;
        Covg sum_covg_forward, sum_covg_reverse;
        boolean cycle_detected;
        
        
        char curr_node_status;
        
        sum_covg_forward = walk_through_graph(db_graph, curr_node, forward, 0, &cycle_detected, &contig_length_forward);
//        printf("length = %d, avg_covg = %f\n",contig_length_forward, avg_covg_forward);
        curr_node_status = curr_node->status;
        
        db_node_set_status(curr_node,none);
        printf("finish forward\n");

        sum_covg_reverse = walk_through_graph(db_graph, curr_node, reverse, 0, &cycle_detected, &contig_length_reverse);
//        printf("length = %d, avg_covg = %f\n",contig_length_reverse, avg_covg_reverse);
        
        printf("finish reverse\n");
        
        if (curr_node_status == special_visited)
            db_node_set_status(curr_node,special_visited);
        else if (db_node_check_status(curr_node,none))
            db_node_set_status(curr_node,visited);
            
            
        int contig_length = contig_length_forward + contig_length_reverse - 1;
        double avg_covg = (sum_covg_reverse + sum_covg_forward - db_node_get_coverage(curr_node, 0))/(double)(contig_length);
        
        printf("length = %d, avg_covg = %f\n",contig_length, avg_covg);
        
        if (contig_length<kmer_size)
            num_small_contig++;
        
        if (num_small_contig > 10)
            break;
        
        
        //set the visited nodes status to keep
        for(i=0;i < db_graph->number_buckets * db_graph->bucket_size;i++)
        {
          if (!db_node_check_for_flag_ALL_OFF(&db_graph->table[i]) && !db_node_check_status(&db_graph->table[i],none))
          {
            if ((contig_length<kmer_size) && (db_node_check_status(&db_graph->table[i],visited) || db_node_check_status(&db_graph->table[i],special_visited)))
                db_node_set_status(&db_graph->table[i],pruned);
            else if (db_node_check_status(&db_graph->table[i],visited))
                db_node_set_status(&db_graph->table[i],keep);
            else if (db_node_check_status(&db_graph->table[i],special_visited))
                db_node_set_status(&db_graph->table[i],special_keep);
            
            
          }  
        }
        
    }
    
    //create a small dB graph to store the possible genome path
    dBGraph * db_graph_small = NULL;
    db_graph_small = hash_table_new(hash_key_bits-5,bucket_size, max_retries, kmer_size);
    if (db_graph_small==NULL)
    {
        die("Giving up - unable to allocate memory for a small hash table\n");
    }
    
    dBNode *next_node, *current_node_in_small_graph;
    Orientation next_orientation;
    Nucleotide nuc, reverse_nuc;
    Edges edge;
    int n;
    BinaryKmer tmp_kmer;
    
    for(i=0;i < db_graph->number_buckets * db_graph->bucket_size;i++)
    {
        if (!db_node_check_for_flag_ALL_OFF(&db_graph->table[i]) && (db_node_check_status(&db_graph->table[i],keep) || db_node_check_status(&db_graph->table[i],special_keep)))
        {
            
            current_node_in_small_graph = hash_table_insert(element_get_key(element_get_kmer(&db_graph->table[i]),db_graph_small->kmer_size, &tmp_kmer),db_graph_small);
            db_node_set_status(current_node_in_small_graph, none);
            
            db_node_update_coverage(current_node_in_small_graph, 0, db_node_get_coverage(&db_graph->table[i], 0));

            //forward
            edge = db_node_get_edge_in_specific_person_or_population(&db_graph->table[i], forward, 0);
            
            for(n=0;n<4;n++)
            {
                if ((edge & 1) == 1)
                {
                    nuc = n;
                  
                    next_node = db_graph_get_next_node_for_specific_person_or_pop(&db_graph->table[i], forward, &next_orientation, nuc, &reverse_nuc, db_graph, 0);
                    if (next_node == NULL)
                        continue;
                    
                    if (db_node_check_status(next_node,keep) || db_node_check_status(next_node,special_keep))
                        db_node_add_labeled_edge(current_node_in_small_graph, forward, nuc, 0);
                }
                
                edge >>= 1;    
            }

            // reverse       
            edge = db_node_get_edge_in_specific_person_or_population(&db_graph->table[i], reverse, 0);
            
            for(n=0;n<4;n++)
            {
                if ((edge & 1) == 1)
                {
                    nuc = n;
                  
                    next_node = db_graph_get_next_node_for_specific_person_or_pop(&db_graph->table[i], reverse, &next_orientation, nuc, &reverse_nuc, db_graph, 0);
                    if (next_node == NULL)
                        continue;
                    
                    if (db_node_check_status(next_node,keep) || db_node_check_status(next_node,special_keep))
                        db_node_add_labeled_edge(current_node_in_small_graph, reverse, nuc, 0);
                }
                
                edge >>= 1;    
            }        
      
        }
    }
    
    db_graph_health_check(true, db_graph_small);
    //dump the small db_graph
//    if (cmd_line->dump_binary==true)
//        db_graph_dump_single_colour_binary_of_colour0(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,
//							db_graph_small, db_graph_info, BINVERSION);
    
    ///////
    
    int count_forward, count_reverse;
    int next_node_in, next_node_out;

    // remove bubble
    for(i=0;i < db_graph_small->number_buckets * db_graph_small->bucket_size;i++)
    {
        if (!db_node_check_for_flag_ALL_OFF(&db_graph_small->table[i]) && db_node_check_status(&db_graph_small->table[i], none))
        {
            
            count_forward = count_num_of_edges_all_colours(&db_graph_small->table[i],forward);
            count_reverse = count_num_of_edges_all_colours(&db_graph_small->table[i],reverse);
            //remove bubble
            
            dBNode* bubble_end_node[2];
            Orientation bubble_end_orient[2];
            double bubble_avg_coverage[2];
            int bubble_n[2];
            
            int bubble_index = 0;
            
            dBNode* start_node;
            Orientation start_orient;
            
            int n_temp;
            Nucleotide nuc_temp, reverse_nuc_temp;
            Edges edge_temp;
            
            
                        
            if (count_forward == 2)
            {
                start_node = &db_graph_small->table[i];
                start_orient = forward;

                
                edge = db_node_get_edge_in_specific_person_or_population(start_node, start_orient, 0);
                    
                for(n=0;n<4;n++)
                {
                    if ((edge & 1) == 1)
                    {
                        int bubble_len = 0;
                        Covg bubble_covg = 0;


                        nuc = n;

                        next_node = db_graph_get_next_node_for_specific_person_or_pop(start_node, start_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);

                        if (next_node == NULL)
                            die("somthing wrong in small db_graph\n");
                        
                        curr_node = next_node;
                        curr_orient = next_orientation;
                            
                        do
                        {
                            db_node_set_status(curr_node, visited);
                            bubble_covg += db_node_get_coverage(curr_node, 0);
                            bubble_len++;
                            
                            edge_temp = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                            for(n_temp=0;n_temp<4;n_temp++)
                            {
                                if ((edge_temp & 1) == 1)
                                {
                                    nuc_temp = n_temp;
                                  
                                    next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc_temp, &reverse_nuc_temp, db_graph_small, 0);
                                    if (next_node == NULL)
                                        continue;
                                    
                                    
                                }
                                
                                edge_temp >>= 1;    
                            }
                            
                            curr_node = next_node;
                            curr_orient = next_orientation;
                            
                            
                        }while(count_num_of_edges_all_colours(curr_node,opposite_orientation(curr_orient)) == 1 && count_num_of_edges_all_colours(curr_node,curr_orient) == 1);
                       
                    
                        
                        bubble_end_node[bubble_index] = curr_node;
                        bubble_end_orient[bubble_index] = curr_orient;
                        bubble_avg_coverage[bubble_index] = (float)bubble_covg/bubble_len;
                        bubble_n[bubble_index] = n;
                        bubble_index++;
                    }
                    
                    edge >>= 1;    
                }
                    
                
                //set smaller bubble to be pruned
                if ((bubble_end_node[0] == bubble_end_node[1]) && (bubble_end_orient[0] == bubble_end_orient[1]))
                {
                    printf("Bubble found! Keep the one with larger coverage!\n");
                    
                    if (bubble_avg_coverage[0] > bubble_avg_coverage[1])
                        nuc = bubble_n[1];
                    else
                        nuc = bubble_n[0];
                        
                    
                    next_node = db_graph_get_next_node_for_specific_person_or_pop(start_node, start_orient, &next_orientation, nuc, &reverse_nuc, db_graph_small, 0);
                    
                    curr_node = next_node;
                    curr_orient = next_orientation;
                    
                    //reset edges
                    db_node_reset_edge(start_node,start_orient, nuc, 0);
                    db_node_reset_edge(next_node,opposite_orientation(next_orientation), reverse_nuc, 0);
                    
                    
                    do
                    {
                        db_node_set_status(curr_node, pruned);
                        edge_temp = db_node_get_edge_in_specific_person_or_population(curr_node, curr_orient, 0);
                        for(n_temp=0;n_temp<4;n_temp++)
                        {
                            if ((edge_temp & 1) == 1)
                            {
                                nuc_temp = n_temp;
                              
                                next_node = db_graph_get_next_node_for_specific_person_or_pop(curr_node, curr_orient, &next_orientation, nuc_temp, &reverse_nuc_temp, db_graph_small, 0);
                                if (next_node == NULL)
                                    continue;
                                
                                
                            }
                            
                            edge_temp >>= 1;    
                        }
                        
                        
                        next_node_in = count_num_of_edges_all_colours(next_node,opposite_orientation(next_orientation));
                        next_node_out = count_num_of_edges_all_colours(next_node,next_orientation);
                        
                        db_node_reset_edge(curr_node,curr_orient, nuc_temp, 0);
                        db_node_reset_edge(next_node,opposite_orientation(next_orientation), reverse_nuc_temp, 0);
                       
                        curr_node = next_node;
                        curr_orient = next_orientation;
                        
                        
                    }while(next_node_in == 1 && next_node_out == 1);
                
                 
                    db_node_set_status(curr_node, visited);
                    db_node_set_status(start_node, visited);
                    
                }
              
            }
//            printf("fwd = %d, rev = %d\n", count_forward, count_reverse);
        }
    }
    
    //health check
    db_graph_health_check(true, db_graph_small);
    hash_table_traverse(&db_node_action_set_status_of_unpruned_to_none, db_graph_small);
    
    //dump the small db_graph
//    if (cmd_line->dump_binary==true)
//        db_graph_dump_single_colour_binary_of_colour0(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,
//							db_graph_small, db_graph_info, BINVERSION);
//    

  // remove small tips
    
    remove_shorter_tips(db_graph_small, reverse);
    remove_shorter_tips(db_graph_small, forward);
    

//    db_graph_health_check(true, db_graph_small);
    
    hash_table_traverse(&db_node_action_set_status_of_unpruned_to_none, db_graph_small);
    //dump the small db_graph
    if (cmd_line->dump_binary==true)
        db_graph_dump_single_colour_binary_of_colour0(cmd_line->output_binary_filename, &db_node_check_status_not_pruned,
							db_graph_small, db_graph_info, BINVERSION);    

}