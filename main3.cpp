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
#include <pbc/pbc.h>

using namespace std::chrono;
using namespace std;
using namespace seal;

int main(int argc, char *argv[]) {

    uint64_t number_of_items = 1 << 16;
    uint64_t size_per_item = 20; // in bytes
    uint32_t N = 2048;

    size_t item_of_layer = 1 << 16;
    size_t number_of_layer = 1;


    // Recommended values: (logt, d) = (12, 2) or (8, 1). 
    uint32_t logt = 12; 
    uint32_t d = 2;

    EncryptionParameters params(scheme_type::BFV);
    PirParams pir_params; 

    // Generates all parameters
    cout << "Main: Generating all parameters" << endl;
    gen_params(number_of_items, size_per_item, N, logt, d, params, pir_params);

    cout << "Main: Initializing the database (this may take some time) ..." << endl;

    // Create test database
    auto db(make_unique<uint8_t[]>(number_of_items * size_per_item));

    // Copy of the database. We use this at the end to make sure we retrieved
    // the correct element.
    auto db_copy(make_unique<uint8_t[]>(number_of_items * size_per_item));

    auto time_generate_s = high_resolution_clock::now();

    mpz_t * mpz_db;
    gmp_randstate_t state;
    mpz_t order;
    gmp_randinit_default(state);
    mpz_init(order);    
    mpz_set_str(order,"730750818665451621361119245571504901405976559617",10);

    mpz_db = (mpz_t *)malloc(sizeof(mpz_t) * number_of_items); 
    size_t byte_count;
    for (size_t i = 0; i < number_of_items; i ++){
        mpz_init(mpz_db[i]);
        mpz_urandomm(mpz_db[i], state, order);
        byte_count = (mpz_sizeinbase (mpz_db[i], 2) + 7 ) / 8;
        mpz_export(db.get() + i * size_per_item + (size_per_item - byte_count),\
                     NULL, 1, sizeof(db.get()[0]), 0, 0, mpz_db[i]);
        mpz_export(db_copy.get() + i * size_per_item + (size_per_item - byte_count),\
                     NULL, 1, sizeof(db.get()[0]), 0, 0, mpz_db[i]);
    }
    
    random_device rd;
    // for (uint64_t i = 0; i < number_of_items; i++) {
    //     for (uint64_t j = 0; j < size_per_item; j++) {
    //         auto val = rd() % 256;
    //         db.get()[(i * size_per_item) + j] = val;
    //         db_copy.get()[(i * size_per_item) + j] = val;
    //     }
    // }

    auto time_generate_e = high_resolution_clock::now();
    auto time_generate_us = duration_cast<microseconds>(time_generate_e - time_generate_s).count();

    // mpz_t mpz_db[number_of_items];
    // for (size_t i = 0; i < number_of_items; i++){
    //     mpz_init(mpz_db[i]);
    //     mpz_import(mpz_db[i], size_per_item, 1, sizeof(uint8_t), 0, 0, db.get() + i * size_per_item);
    // }

    auto time_vc_s = high_resolution_clock::now();

    pairing_t pairing;
    FILE* fp;
    fp = fopen("/home/liwenming/桌面/pbc-0.5.14/param/a.param","r");
    char param[1024];
    size_t count = fread(param, 1, 1024, fp);

    pairing_init_set_buf(pairing, param, count);

    element_t g1, g2, gt;

    element_init_G1(g1, pairing);
    element_init_G2(g2, pairing);
    element_init_GT(gt, pairing);

    element_random(g1);
    element_random(g2);
    element_pairing(gt, g1 ,g2);

    element_t* g1_t;
    element_t* g1_tt;
    element_t* g2_t;
    
    g1_t = (element_t *)malloc(sizeof(element_t) * item_of_layer);
    g1_tt = (element_t *)malloc(sizeof(element_t) * (item_of_layer - 1));
    g2_t = (element_t *)malloc(sizeof(element_t) * item_of_layer);

    mpz_t alpha;
    mpz_init(alpha);


    mpz_urandomm(alpha, state, order);

    for (size_t i = 0; i < item_of_layer; i ++){
        element_init_G1(g1_t[i],pairing);
        element_init_G2(g2_t[i],pairing);
    }
    for (size_t i = 0; i < item_of_layer - 1; i ++){
        element_init_G1(g1_tt[i],pairing);
    }

    element_pp_t g1_pp,g2_pp,gt_pp;
    element_pp_init(g1_pp,g1);
    element_pp_init(g2_pp,g2);
    element_pp_init(gt_pp,gt);
    
    mpz_t* a;
    a = (mpz_t *)malloc(sizeof(mpz_t) * item_of_layer);
    mpz_init(a[0]);
    mpz_set(a[0],alpha);
    for (size_t i = 1; i < item_of_layer; i ++ ){
        mpz_init(a[i]);
        mpz_mul(a[i],a[i - 1],alpha);
        mpz_mod(a[i],a[i],order);
    }
    for (size_t i = 0; i < item_of_layer; i ++){
        element_pp_pow(g1_t[i],a[i],g1_pp);
        element_pp_pow(g2_t[i],a[i],g2_pp);
    }

    mpz_t temp;
    mpz_init(temp);
    for (size_t i = 0; i < item_of_layer - 1; i ++){
        mpz_mul(temp, a[item_of_layer - 1], a[i + 1]);
        mpz_mod(temp, temp, order);
        element_pp_pow(g1_tt[i], temp, g1_pp);        
    }
    
    mpz_clear(alpha);

    auto time_vc_m = high_resolution_clock::now();
    auto time_vc_us = duration_cast<microseconds>(time_vc_m - time_vc_s).count();
/* Commit(m) */

    // element_t ** C;
    // C = (element_t **)malloc(sizeof(element_t *) * number_of_layer);
    // mpz_t * target = mpz_db;

    // for (size_t i = 0; i < number_of_layer ; i ++){
    //     size_t layer_size = number_of_items / (size_t)pow(item_of_layer, (i + 1));
    //     C[i] = (element_t *)malloc(sizeof(element_t) * layer_size);
    //     for (size_t j = 0; j <  layer_size; j ++){
    //         element_init_G1(C[i][j], pairing);
    //         mpz_t sum;
    //         mpz_init(sum);
    //         mpz_set_str(sum,"0",10);
    //         for (){

    //         }
    //         mpz_clear(sum);
    //     }
    // }

    element_t C;
    mpz_t sum;
    mpz_init(sum);
    mpz_set_str(sum,"0",10);
    element_init_G1(C,pairing);
    
    for (size_t i = 0; i < item_of_layer; i ++){
        mpz_mul(temp, mpz_db[i], a[i]);
        mpz_mod(temp, temp, order);
        mpz_add(sum, sum, temp);
        mpz_mod(sum, sum, order);
    }
    element_pp_pow(C,sum, g1_pp);

    auto time_vc_e = high_resolution_clock::now();
    auto time_vc2_us = duration_cast<microseconds>(time_vc_e - time_vc_m).count();


    /* prove */

    element_t * proof;
    proof = (element_t *)malloc(sizeof(element_t) * item_of_layer);
    for (size_t i = 0; i < item_of_layer; i ++){
        mpz_set_str(sum, "0", 10);
        element_init_G1(proof[i],pairing);
        for (size_t j = 0; j < item_of_layer; j ++){
            if (j != i){
                mpz_mul(temp, mpz_db[j], a[j]);
                mpz_mod(temp, temp, order);
                mpz_add(sum, sum, temp);
                mpz_mod(sum, sum, order);
            }
        }
        mpz_mul(sum, sum, a[item_of_layer - i - 1]);
        element_pp_pow(proof[i], sum, g1_pp);
    }

    auto time_vc_ee = high_resolution_clock::now();
    auto time_vc3_us = duration_cast<microseconds>(time_vc_ee - time_vc_e).count();



    // Initialize PIR Server
    cout << "Main: Initializing server and client" << endl;
    PIRServer server(params, pir_params);

    // Initialize PIR client....
    PIRClient client(params, pir_params);
    GaloisKeys galois_keys = client.generate_galois_keys();

    // Set galois key for client with id 0
    cout << "Main: Setting Galois keys...";
    server.set_galois_key(0, galois_keys);

    // Measure database setup
    auto time_pre_s = high_resolution_clock::now();
    server.set_database(move(db), number_of_items, size_per_item);
    server.preprocess_database();
    cout << "Main: database pre processed " << endl;
    auto time_pre_e = high_resolution_clock::now();
    auto time_pre_us = duration_cast<microseconds>(time_pre_e - time_pre_s).count();

    // Choose an index of an element in the DB
    uint64_t ele_index = rd() % number_of_items; // element in DB at random position
    uint64_t index = client.get_fv_index(ele_index, size_per_item);   // index of FV plaintext
    uint64_t offset = client.get_fv_offset(ele_index, size_per_item); // offset in FV plaintext
    cout << "Main: element index = " << ele_index << " from [0, " << number_of_items -1 << "]" << endl;
    cout << "Main: FV index = " << index << ", FV offset = " << offset << endl; 

    // Measure query generation
    auto time_query_s = high_resolution_clock::now();
    PirQuery query = client.generate_query(index);
    auto time_query_e = high_resolution_clock::now();
    auto time_query_us = duration_cast<microseconds>(time_query_e - time_query_s).count();
    cout << "Main: query generated" << endl;

    //To marshall query to send over the network, you can use serialize/deserialize:
    //std::string query_ser = serialize_query(query);
    //PirQuery query2 = deserialize_query(d, 1, query_ser, CIPHER_SIZE);

    // Measure query processing (including expansion)
    auto time_server_s = high_resolution_clock::now();
    PirReply reply = server.generate_reply(query, 0);
    auto time_server_e = high_resolution_clock::now();
    auto time_server_us = duration_cast<microseconds>(time_server_e - time_server_s).count();

    // Measure response extraction
    auto time_decode_s = chrono::high_resolution_clock::now();
    Plaintext result = client.decode_reply(reply);
    auto time_decode_e = chrono::high_resolution_clock::now();
    auto time_decode_us = duration_cast<microseconds>(time_decode_e - time_decode_s).count();

    // Convert from FV plaintext (polynomial) to database element at the client
    vector<uint8_t> elems(N * logt / 8);
    coeffs_to_bytes(logt, result, elems.data(), (N * logt) / 8);
    // Check that we retrieved the correct element
    for (uint32_t i = 0; i < size_per_item; i++) {
        if (elems[(offset * size_per_item) + i] != db_copy.get()[(ele_index * size_per_item) + i]) {
            cout << "Main: elems " << (int)elems[(offset * size_per_item) + i] << ", db "
                 << (int) db_copy.get()[(ele_index * size_per_item) + i] << endl;
            cout << "Main: PIR result wrong!" << endl;
            return -1;
        }
    }

    // Output results
    cout << "Main: PIR result correct!" << endl;

    cout << "Main: Database generate time: " << time_generate_us / 1000 << " ms" << endl;
    cout << "Main: Vector commitment pre time: " << time_vc_us / 1000 << " ms" << endl;
    cout << "Main: Vector commitment commit time: " << time_vc2_us / 1000 << " ms" << endl;
    cout << "Main: Vector commitment prove(simple) time: " << time_vc3_us / 1000 << " ms" << endl;
    
    cout << "Main: PIRServer pre-processing time: " << time_pre_us / 1000 << " ms" << endl;
    cout << "Main: PIRClient query generation time: " << time_query_us / 1000 << " ms" << endl;
    cout << "Main: PIRServer reply generation time: " << time_server_us / 1000 << " ms"
         << endl;
    cout << "Main: PIRClient answer decode time: " << time_decode_us / 1000 << " ms" << endl;
    cout << "Main: Reply num ciphertexts: " << reply.size() << endl;
    
    return 0;
}
