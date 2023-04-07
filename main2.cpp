#include "pir.hpp"
#include "pir_client.hpp"
#include "pir_server.hpp"
#include <seal/seal.h>
#include <chrono>
#include <memory>
#include <random>
#include <cstdint>
#include <cstddef>
#include <math.h>
#include <gmp.h>
#include <memory>
#include <pbc/pbc.h>

using namespace std::chrono;
using namespace std;
using namespace seal;

int main(int argc, char *argv[]) {

    uint64_t number_of_items = 1 << 20;
    uint32_t N = 2048;
    size_t size_per_item = 20;

    size_t number_of_layer = 3;
    size_t item_of_commit = 1 << 7;

    // Recommended values: (logt, d) = (12, 2) or (8, 1). 
    uint32_t logt = 12; 
    uint32_t d = 2;

    FILE* fp;
    fp = fopen("/home/liwenming/桌面/pbc-0.5.14/param/a.param","r");
    pairing_t pairing;
    char param[1024];
    size_t count = fread(param, 1, 1024, fp);
    pairing_init_set_buf(pairing, param, count);

    mpz_t order;
    mpz_init(order);    
    mpz_set_str(order,"730750818665451621361119245571504901405976559617",10);


//_____________________________________________________________________

    /* generate the random database */
    mpz_t ** mpz_db;
    unique_ptr<uint8_t[]>  db[number_of_layer];
    
    gmp_randstate_t state;
    gmp_randinit_default(state);

    mpz_db = (mpz_t **)malloc(sizeof(mpz_t *) * number_of_layer);
    mpz_db[0] = (mpz_t *)malloc(sizeof(mpz_t) * number_of_items); 
    db[0] = make_unique<uint8_t[]>(number_of_items * size_per_item);   
    unique_ptr<uint8_t[]> db_copy = make_unique<uint8_t[]>(number_of_items * size_per_item);
    size_t byte_count;
    for (size_t i = 0; i < number_of_items; i ++){
        mpz_init(mpz_db[0][i]);
        mpz_urandomm(mpz_db[0][i], state, order);
        byte_count = (mpz_sizeinbase (mpz_db[0][i], 2) + 7 ) / 8;
        mpz_export(db[0].get() + i * size_per_item + (size_per_item - byte_count),\
                     NULL, 1, sizeof(db[0].get()[0]), 0, 0, mpz_db[0][i]);
        mpz_export(db_copy.get() + i * size_per_item + (size_per_item - byte_count),\
                     NULL, 1, sizeof(db[0].get()[0]), 0, 0, mpz_db[0][i]);
    }

cout << "size_per_item = " << size_per_item << "\n";
cout << "mpz_db[" << 0 << "] has size " << number_of_items << "\n";
cout << "db[" << 0 << "] has size " << number_of_items  * size_per_item << "\n";



   /* initialize the pairing, based on item_of_commit */
cout << "___start initialize____\n";
auto time_init_s = high_resolution_clock::now();
    element_t g1, g2, gt;

    element_init_G1(g1, pairing);
    element_init_G2(g2, pairing);
    element_init_GT(gt, pairing);

    element_random(g1);
    element_random(g2);
    element_pairing(gt, g1 ,g2);

    size_t size_per_element_item = element_length_in_bytes(g1);
cout << "size_per_element_item = " << size_per_element_item << "\n";

    element_t* g1_t;
    element_t* g1_tt;
    element_t* g2_t;
    
     // g1_t[i] = g1 ^ (alpha^i); same as g2_t[i]; 
    // keep g1_t[0],g1_t[N+1],g2_t[0] should not reveal
    // todo: now is reveal
    g1_t = (element_t *)malloc(sizeof(element_t) * (2 * item_of_commit + 1) );
    g2_t = (element_t *)malloc(sizeof(element_t) * (item_of_commit + 1));

    mpz_t alpha;
    mpz_init(alpha);
    mpz_urandomm(alpha, state, order);

    for (size_t i = 1; i <= 2 * item_of_commit; i ++){
        element_init_G1(g1_t[i],pairing);
    }
    for (size_t i = 1; i <= item_of_commit; i ++){
        element_init_G2(g2_t[i],pairing);
    }

    element_pp_t g1_pp,g2_pp,gt_pp;
    element_pp_init(g1_pp,g1);
    element_pp_init(g2_pp,g2);
    element_pp_init(gt_pp,gt);

    // a[i] = alpha^i ; 
    // a[0] is useless, a can not reveal;
    // todo: a is now reveal  
    mpz_t* a;
    a = (mpz_t *)malloc(sizeof(mpz_t) * (item_of_commit + 1));
    mpz_init(a[0]);
    mpz_set_str(a[0],"1",10);
    for (size_t i = 1; i <= item_of_commit; i ++ ){
        mpz_init(a[i]);
        mpz_mul(a[i],a[i - 1],alpha);
        mpz_mod(a[i],a[i],order);
    }
    mpz_t an;
    mpz_init(an);
    for (size_t i = 1; i <= item_of_commit; i ++){
        mpz_mul(an,a[i],a[number_of_layer]);
        mpz_mod(an,an,order);
        element_pp_pow(g1_t[i],a[i],g1_pp);
        if (i!=1){element_pp_pow(g1_t[i + item_of_commit], an, g1_pp);}
        element_pp_pow(g2_t[i],a[i],g2_pp);
    }
    
    mpz_clear(alpha);
// _________________________________________________
 

/* Commit(m) */

    std::cout << "____________start commit___________\n";
    auto time_commit_s = high_resolution_clock::now();

    element_t ** C;
    C = (element_t **)malloc(sizeof(element_t *) * (number_of_layer + 1));

    size_t item_of_layer[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++){
        item_of_layer[i] = number_of_items / (size_t)pow(item_of_commit, i);
    }
    size_t item_of_final_commit = item_of_layer[number_of_layer - 1];

    mpz_t * target_db;
    mpz_t sum;
    mpz_t mul;
    mpz_t residual;
    mpz_init(sum);
    mpz_init(mul);
    mpz_init(residual);
            
    for (size_t i = 1; i < number_of_layer ; i ++){
        C[i] = (element_t *)malloc(sizeof(element_t) * item_of_layer[i]);
        mpz_db[i] = (mpz_t *)malloc(sizeof(mpz_t) * item_of_layer[i]);
        db[i] = make_unique<uint8_t[]>(item_of_layer[i] * size_per_element_item);

cout << "mpz_db[" << i << "] has size " << item_of_layer[i] << "\n";
cout << "db[" << i << "] has size " << item_of_layer[i] * size_per_element_item << "\n";

        for (size_t j = 0; j < item_of_layer[i]; j ++){
            //  cout << "j = " << j << "\n";
            // cout << "here is j = " << j << ",";
            target_db = mpz_db[i - 1] + j * item_of_commit;
            element_init_G1(C[i][j], pairing);
// cout << "j = " << j << "\n";                 
            mpz_set_str(sum,"0",10);
            // C = g1 ^ (mTa)
            for (size_t k = 1; k <= item_of_commit; k ++){
                mpz_mod(residual, target_db[k - 1], order);
                mpz_mul(mul,residual, a[k]);
                mpz_mod(mul, mul, order);
                mpz_add(sum, sum, mul);
                mpz_mod(sum, sum, order);
            }

            element_pp_pow(C[i][j], sum ,g1_pp); 
            element_to_bytes(db[i].get() + j * size_per_element_item, C[i][j]);
            mpz_import(mpz_db[i][j], size_per_element_item, 1, sizeof(db[i].get()[0]),\
                         0, 0, db[i].get() + j * size_per_element_item);            
        }
    }
    // generate the public value: 
cout << "se\n";
    C[number_of_layer] = (element_t *)malloc(sizeof(element_t));
    target_db = mpz_db[number_of_layer - 1];
    element_init_G1(C[number_of_layer][0],pairing);
    
    mpz_set_str(sum,"0",10);
    for (size_t k = 1; k <= item_of_final_commit; k ++){
        mpz_mod(residual, target_db[k - 1], order);
        mpz_mul(mul,residual, a[k]);
        mpz_mod(mul, mul, order);
        mpz_add(sum, sum, mul);
        mpz_mod(sum, sum, order);    
    }
    element_pp_pow(C[number_of_layer][0],sum, g1_pp);

    auto time_commit_e = high_resolution_clock::now();
    auto time_commit_us = duration_cast<microseconds>(time_commit_e - time_commit_s).count();
    std::cout << "_______end commit_____________\n";
    std::cout << "commit take time " << time_commit_us / 1000 << " ms\n";


    /* prove */


    auto time_prove_s = high_resolution_clock::now();
    std::cout << "_______start prove_____________\n";
    mpz_t sub;
    element_t ** proof;
    unique_ptr<uint8_t[]> proof_db[number_of_layer];
    mpz_init(sub);

    proof = (element_t **)malloc(sizeof(element_t*) * number_of_layer);
    for (size_t i = 0; i < number_of_layer; i ++){
        proof_db[i] = make_unique<uint8_t[]>(item_of_layer[i] * size_per_element_item);
        size_t commit_of_this_layer = (item_of_layer[i] + item_of_commit - 1) / item_of_commit;
cout << "i = " << i << " , commit of this layer = " << commit_of_this_layer << "\n";
        proof[i] = (element_t *)malloc(sizeof(element_t) * item_of_layer[i]);
        for (size_t j = 0; j < commit_of_this_layer; j ++){
            if (i != number_of_layer - 1){
                mpz_set_str(sum,"0",10); 
                for (size_t l = 0; l < item_of_commit; l ++){
                    mpz_mod(residual,mpz_db[i][j * item_of_commit + l], order);
                    mpz_mul(mul, residual, a[l + 1]);
                    mpz_mod(mul, mul, order);
                    mpz_add(sum, sum, mul);
                    mpz_mod(sum, sum, order);
                }


                for (size_t k = 0; k < item_of_commit; k ++){
                    element_init_G1(proof[i][j * item_of_commit + k],pairing);
                    // proof_i  = g_1^(a^(N+1-i) * m[-i]Ta[-i])
                    // todo: here use a to generate proof, server can also use g1 to generate  
                    mpz_mod(residual,mpz_db[i][j * item_of_commit + k], order);
                    mpz_mul(mul, residual, a[k + 1]);
                    mpz_mod(mul, mul, order);
                    mpz_sub(sub, sum, mul);

                    mpz_mul(mul, sub, a[item_of_commit - k]);
                    mpz_mod(mul, mul, order);

                    auto pairing_s = high_resolution_clock::now();
                    element_pp_pow(proof[i][j * item_of_commit + k], mul ,g1_pp);
                    auto pairing_e = high_resolution_clock::now();
                    
                }
            }
            else{
                cout << "final layer\n";
                mpz_set_str(sum,"0",10);
                for (size_t l = 0; l < item_of_final_commit; l ++){
                    mpz_mod(residual,mpz_db[i][j * item_of_final_commit + l], order);
                    mpz_mul(mul, residual, a[l + 1]);
                    mpz_mod(mul, mul, order);
                    mpz_add(sum, sum, mul);
                    mpz_mod(sum, sum, order);
                }



                for (size_t k = 0; k < item_of_final_commit; k ++){
                    element_init_G1(proof[i][j * item_of_final_commit + k],pairing);
                    // proof_i  = g_1^(a^(N+1-i) * m[-i]Ta[-i])
                    // todo: here use a to generate proof, server can also use g1 to generate  
                    mpz_mod(residual,mpz_db[i][j * item_of_final_commit + k], order);
                    mpz_mul(mul, residual, a[k + 1]);
                    mpz_mod(mul, mul, order);
                    mpz_sub(sub, sum, mul);

                    mpz_mul(mul, sub, a[item_of_final_commit - k]);
                    mpz_mod(mul, mul, order);
                    element_pp_pow(proof[i][j * item_of_final_commit + k], mul ,g1_pp);
                }
            }
        }
    }
    for (size_t i = 0; i < number_of_layer; i ++){
        for (size_t j = 0; j < item_of_layer[i]; j ++){
            element_to_bytes(proof_db[i].get() + j * size_per_element_item, proof[i][j]);
        }
    }
    auto time_prove_e = high_resolution_clock::now();
    auto time_prove_us = duration_cast<microseconds>(time_prove_e - time_prove_s).count();
    std::cout << "_______end prove___________\n";
    std::cout << "prove take time " << time_prove_us / 1000 << " ms\n";


//——————————————————————————————————————————————————————

    uint64_t ele_index = 132;
// ———————————————————————————————————————————————————————
    mpz_t client_db_query[number_of_layer];
    element_t client_proof_query[number_of_layer];
    element_t client_commit_query[number_of_layer];
    size_t ele_index_of_layer[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++){
        mpz_init(client_db_query[i]);
        element_init_G1(client_proof_query[i], pairing);
        element_init_G1(client_commit_query[i], pairing);
        ele_index_of_layer[i] = ele_index / (size_t)pow(item_of_commit,i);
        cout << "query index of layer " << i << " is " << ele_index_of_layer[i]  <<"\n";
    }

    EncryptionParameters params(scheme_type::BFV);
    PirParams pir_params;
    allocator<PIRServer> alloc_server;
    allocator<PIRClient> alloc_client;
    PIRServer *server = alloc_server.allocate(number_of_layer);
    PIRClient *client = alloc_client.allocate(number_of_layer);

    allocator<PIRServer> alloc_server_proof;
    allocator<PIRClient> alloc_client_proof;
    PIRServer *server_proof = alloc_server_proof.allocate(number_of_layer);
    PIRClient *client_proof = alloc_client_proof.allocate(number_of_layer);
    
    for (size_t i = 0; i < number_of_layer; i ++){
        size_t size;
        if (i == 0){size = size_per_item;} else {size = size_per_element_item;} 
        gen_params(item_of_layer[i], size, N, logt, d, params, pir_params);
        
        alloc_server.construct(server + i, params, pir_params);
        alloc_client.construct(client + i, params, pir_params);
        GaloisKeys galois_keys = client[i].generate_galois_keys();
        server[i].set_galois_key(0, galois_keys);

        server[i].set_database(move(db[i]), item_of_layer[i], size);
        server[i].preprocess_database();
    }
    for (size_t i = 0; i < number_of_layer; i ++){
        gen_params(item_of_layer[i], size_per_element_item, N ,logt, d, params, pir_params);

        alloc_server_proof.construct(server_proof + i, params, pir_params);
        alloc_client_proof.construct(client_proof + i, params, pir_params);
        GaloisKeys galois_keys = client_proof[i].generate_galois_keys();
        server_proof[i].set_galois_key(0, galois_keys);

        server_proof[i].set_database(move(proof_db[i]), item_of_layer[i], size_per_element_item);
        server_proof[i].preprocess_database();
    }

    auto time_init_e = high_resolution_clock::now();
    auto time_init_us = duration_cast<microseconds>(time_init_e - time_init_s).count();
    cout << "____________________\n";
    std::cout << "Initialization total take time " << time_init_us / 1000 << " ms\n";
    cout << "____________________\n";

// Query
    auto time_query_s = high_resolution_clock::now();
    PirQuery query[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++)
    {
        size_t size;
        if (i == 0){size = size_per_item;} else {size = size_per_element_item;} 
        uint64_t index = client[i].get_fv_index(ele_index_of_layer[i], size);   

        query[i] = client[i].generate_query(index);
    }
    
    PirQuery query_proof[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++)
    {
        uint64_t index = client_proof[i].get_fv_index(ele_index_of_layer[i], size_per_element_item);   

        query_proof[i] = client_proof[i].generate_query(index);
    }
    auto time_query_e = high_resolution_clock::now();
    auto time_query_us = duration_cast<microseconds>(time_query_e - time_query_s).count();
    cout << "____________________\n";
    std::cout << "Query total take time " << time_query_us / 1000 << " ms\n";
    cout << "____________________\n";

// Reply
    auto time_reply_s = high_resolution_clock::now();
    PirReply reply[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++)
    {
        reply[i] = server[i].generate_reply(query[i], 0);
    }
    PirReply reply_proof[number_of_layer];
    for (size_t i = 0; i < number_of_layer; i ++)
    {
        reply_proof[i] = server_proof[i].generate_reply(query_proof[i], 0);
    }
    auto time_reply_e = high_resolution_clock::now();
    auto time_reply_us = duration_cast<microseconds>(time_reply_e - time_reply_s).count();
    cout << "____________________\n";
    std::cout << "Reply total take time " << time_reply_us / 1000 << " ms\n";
    cout << "____________________\n";


// Decode
    auto time_decode_s = high_resolution_clock::now();

    for (size_t i = 0; i < number_of_layer; i ++)
    {
        size_t size;
        if (i == 0){size = size_per_item;} else {size = size_per_element_item;} 
        uint64_t index = client[i].get_fv_index(ele_index_of_layer[i], size);   
        uint64_t offset = client[i].get_fv_offset(ele_index_of_layer[i], size); 

        Plaintext result = client[i].decode_reply(reply[i]);
        uint8_t elems[N * logt / 8];
        coeffs_to_bytes(logt, result, elems, (N * logt) / 8);
        mpz_import(client_db_query[i], size, 1, sizeof(uint8_t), 0, 0, elems + offset * size);

        if (mpz_cmp(client_db_query[i], mpz_db[i][ele_index_of_layer[i]]) == 0){
            cout << "PIR of db[" << i << "] result correct.\n";
        }
        else{
            cout << "PIR of db[" << i << "] result wrong.\n";
        }
        if (i > 0){
            element_from_bytes( client_commit_query[i - 1] , elems + offset * size);
            if ( element_cmp(client_commit_query[i - 1], C[i][ele_index_of_layer[i]]) == 0){
                cout << "PIR of commit[" << i << "] result correct.\n";
            }
            else {
                cout << "PIR of commit[" << i << "] result wrong.\n";
            }
        }
    }
    element_set(client_commit_query[number_of_layer - 1],C[number_of_layer][0]);


    for (size_t i = 0; i < number_of_layer; i ++)
    {
        
        uint64_t index = client_proof[i].get_fv_index(ele_index_of_layer[i], size_per_element_item);   
        uint64_t offset = client_proof[i].get_fv_offset(ele_index_of_layer[i], size_per_element_item); 

        Plaintext result = client_proof[i].decode_reply(reply_proof[i]);
        uint8_t elems[N * logt / 8];
        coeffs_to_bytes(logt, result, elems, (N * logt) / 8);

        element_from_bytes( client_proof_query[i] , elems + offset * size_per_element_item);
        if (element_cmp( client_proof_query[i], proof[i][ele_index_of_layer[i]] ) == 0){
            cout << "PIR of proof[" << i << "] result correct.\n";
        }
        else {
            cout << "PIR of proof[" << i << "] result wrong.\n";
        }
    }
    auto time_decode_e = high_resolution_clock::now();
    auto time_decode_us = duration_cast<microseconds>(time_decode_e - time_decode_s).count();
    cout << "____________________\n";
    std::cout << "Decode total take time " << time_decode_us / 1000 << " ms\n";
    cout << "____________________\n";

     
// verify  e(C,g2^(alpha^( N + 1 - i))) ?= e(proof, g2) * gt^(alpha^(N + 1) * mi)
    auto time_verify_s = high_resolution_clock::now();

    mpz_t aNp1;
    mpz_init(aNp1);
    mpz_mul(aNp1, a[1], a[item_of_commit]);
    mpz_mod(aNp1,aNp1, order);

    bool flag = true;
    element_t left,right;
    element_t right_gt;
    element_init_GT(left, pairing);
    element_init_GT(right, pairing);
    element_init_GT(right_gt, pairing);

    for (size_t i = 0; i < number_of_layer - 1; i ++){
        size_t query_of_commit = ele_index_of_layer[i] % item_of_commit;
        element_pairing(left, client_commit_query[i], g2_t[item_of_commit - query_of_commit]);
        element_pairing(right, client_proof_query[i], g2);
        
        mpz_mul(mul, aNp1, client_db_query[i]);
        mpz_mod(mul, mul, order);
        element_pp_pow(right_gt,mul,gt_pp);
        element_mul(right, right, right_gt);
        if (element_cmp(left, right) != 0){
            std::cout << "error occur at";
            std::cout << " layer "<< i;
            flag = false;
            break;
        }   
    }
    size_t query_of_commit = ele_index_of_layer[number_of_layer - 1] % item_of_final_commit;
    element_pairing(left, client_commit_query[number_of_layer - 1], g2_t[item_of_final_commit - query_of_commit]);
    element_pairing(right, client_proof_query[number_of_layer - 1], g2);
    
    mpz_t af;
    if (item_of_final_commit == item_of_commit){mpz_set(af,aNp1);}else{mpz_set(af,a[item_of_final_commit + 1]);}
    mpz_mul(mul, af, client_db_query[number_of_layer - 1]);
    mpz_mod(mul, mul, order);
    element_pp_pow(right_gt,mul,gt_pp);
    element_mul(right, right, right_gt);
    if (element_cmp(left, right) != 0){
        std::cout << "error occur at";
        std::cout << " layer "<< number_of_layer - 1;
        flag = false;
    }   


    auto time_verify_e = high_resolution_clock::now();
    auto time_verify_us = duration_cast<microseconds>(time_verify_e - time_verify_s).count();
    cout << "____________________\n";
    std::cout << "verify total take time " << time_verify_us / 1000 << " ms\n";
    cout << "____________________\n";


    if (flag){
        cout << "______________________________\n";
        std::cout << "Initialize total take time " << time_init_us / 1000 << " ms\n";
        std::cout << "Query total take time " << time_query_us / 1000 << " ms\n";
        std::cout << "Reply total take time " << time_reply_us / 1000 << " ms\n";
        std::cout << "Decode total take time " << time_decode_us / 1000 << " ms\n";
        std::cout << "Verify total take time " << time_verify_us / 1000 << " ms\n";
        printf("verify result correct\n");
        cout << (int)time_init_us << "\n";
    }
    else{
        printf("false\n");
        element_printf("%B\n",left);
        element_printf("%B\n", right);
        element_cmp(left, right);
    }
   size_t count1 = 0; 
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      for (size_t j = 0; j < query_proof[i].size(); j ++){
         count1 += query[i][j].size();
      }
   }
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      for (size_t j = 0; j < query_proof[i].size(); j ++){
         count1 += query_proof[i][j].size();
      }
   }
   cout << "query size " << count1 * N * logt / 8 << " Bytes\n";
   
   size_t count2 = 0;
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      count2 += reply[i].size();
   }
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      count2 += reply_proof[i].size();
   }
   cout << "reply size " << count2 * N * logt / 8 << " Bytes\n"; 


    // mpz_t aNp1;
    // mpz_init(aNp1);
    // mpz_mul(aNp1, a[1], a[item_of_commit]);
    // mpz_mod(aNp1,aNp1, order);

    // bool flag = true;
    // element_t left,right;
    // element_t right_gt;
    // element_init_GT(left, pairing);
    // element_init_GT(right, pairing);
    // element_init_GT(right_gt, pairing);


    // for (size_t i = 0; i < number_of_layer - 1; i ++ ){
        
    //     size_t item_of_this_layer = number_of_items / (size_t)pow(item_of_commit, i);
    //     size_t commit_of_this_layer = item_of_this_layer / item_of_commit;

    //     for (size_t j = 0; j < commit_of_this_layer; j ++){
    //         // targetC = C[i + 1][j]
    //         for (size_t k = 0; k < item_of_commit; k ++){
    //             element_pairing(left, C[i + 1][j], g2_t[item_of_commit - k]);
    //             element_pairing(right, proof[i][j * item_of_commit + k], g2);

    //             mpz_mul(mul, aNp1, mpz_db[i][j * item_of_commit + k]);
    //             mpz_mod(mul, mul, order);
    //             element_pp_pow(right_gt,mul, gt_pp);
    //             element_mul(right, right, right_gt);
    //             if (element_cmp(left, right) != 0){
    //                 std::cout << "error occur at";
    //                 std::cout << " layer "<< i;
    //                 std::cout <<  ", index = " << j * item_of_commit + k << "\n";
    //                 flag = false;
    //                 break;
    //             }
    //         }
    //         if (!flag){
    //             break;
    //         }
    //     }
    //     if (!flag){
    //         break;
    //     }
    // }




    return 0;
}
