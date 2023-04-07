using namespace std;
#define SHA256_DIGEST_LENGTH 32
#include <openssl/sha.h>
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

using namespace std::chrono;
using namespace std;
using namespace seal;

void sha256(unsigned char * str, size_t size, unsigned char * hash)
{
   // unsigned char * hash = new unsigned char(SHA256_DIGEST_LENGTH);
   SHA256_CTX sha256;
   SHA256_Init(&sha256);
   SHA256_Update(&sha256, str, size);
   SHA256_Final(hash, &sha256);
}


int main() {

    uint64_t number_of_items = 1 << 12;
    size_t number_of_layer = 13;
    
    uint64_t size_per_item = 20; // in bytes
    uint32_t N = 2048;

    // Recommended values: (logt, d) = (12, 2) or (8, 1). 
    uint32_t logt = 12; 
    uint32_t d = 2;
// ________________________________________________________________


    unique_ptr<uint8_t[]> db[number_of_layer];
    unique_ptr<uint8_t[]> db_copy[number_of_layer];

    db[0] = make_unique<uint8_t[]>(number_of_items * size_per_item);
    db_copy[0] = make_unique<uint8_t[]>(number_of_items * size_per_item);
    size_t item_of_layer[number_of_layer];
    item_of_layer[0] = number_of_items;

    for (size_t i = 1; i < number_of_layer; i ++){
       item_of_layer[i] = item_of_layer[i - 1] >> 1;
       db[i] = (make_unique<uint8_t[]>(item_of_layer[i] * SHA256_DIGEST_LENGTH ));
       db_copy[i] = (make_unique<uint8_t[]>(item_of_layer[i] * SHA256_DIGEST_LENGTH));
    }
    random_device rd;
    for (uint64_t i = 0; i < number_of_items; i++) {
        for (uint64_t j = 0; j < size_per_item; j++) {
            auto val = rd() % 256;
            db[0].get()[(i * size_per_item) + j] = val;
            db_copy[0].get()[(i * size_per_item) + j] = val;
        }
    }

   auto time_init_s = high_resolution_clock::now();

   for (size_t i = 1; i < number_of_layer; i ++){
      size_t item_len;
      if (i == 1){item_len = 2 * size_per_item;}else{item_len = 2 * SHA256_DIGEST_LENGTH;}
      
      for (size_t j = 0; j < item_of_layer[i] ; j ++){
         unsigned char data[item_len];
         unsigned char hash[SHA256_DIGEST_LENGTH];
         for (size_t byte_index = 0; byte_index < item_len; byte_index++){
            data[byte_index] = db[i - 1].get()[j * item_len + byte_index];
         }
         sha256(data,item_len,hash);
         if (i == 1 & j == 0){
            for (size_t k = 0; k < item_len; k ++){
               cout << (int)data[k] << " ";
            }
            cout << "\n";
            for (size_t k = 0; k < SHA256_DIGEST_LENGTH; k ++){
               cout << (int)hash[k] << " ";
            }
            cout << "\n";
         }

         for (size_t byte_index = 0; byte_index < SHA256_DIGEST_LENGTH; byte_index ++){
            db[i].get()[j * SHA256_DIGEST_LENGTH + byte_index] = hash[byte_index];
            db_copy[i].get()[j * SHA256_DIGEST_LENGTH + byte_index] = hash[byte_index];
         }
         if (i == 1 & j == 0){
            for (size_t k = 0; k < item_len; k ++){
               cout << (int)db[i].get()[k] << " ";
            }
            cout << "\n";
            for (size_t k = 0; k < SHA256_DIGEST_LENGTH; k ++){
               cout << (int)db_copy[i].get()[k] << " ";
            }
            cout << "\n";
         }        
      }
   }

// output the database db
   
   // for (size_t i = 0; i < number_of_layer; i ++){
   //    size_t size;
   //    cout << " i = " << i <<"\n";
   //    if (i == 0){size = size_per_item;}else{size =  SHA256_DIGEST_LENGTH;}
   //    for (size_t j = 0; j < item_of_layer[i]; j ++){
   //       cout << "j = " << j << ": ";
   //          for (size_t k = 0; k < size; k ++){
   //             cout << (int)db_copy[i].get()[j*size + k] <<" ";
   //          }
   //       cout << "\n";
   //    }
   // }

cout << "here in sha256.cpp \n";



//__________________________________

   uint64_t ele_index = 1;
//____________________________________
   

   EncryptionParameters params(scheme_type::BFV);
   PirParams pir_params;
cout << "here2  in sha256.cpp \n";
 
   // Init the server and client
   allocator<PIRServer> alloc_server;
   allocator<PIRClient> alloc_client;
   PIRServer *server = alloc_server.allocate(number_of_layer);
   PIRClient *client = alloc_client.allocate(number_of_layer);
   size_t ele_index_of_layer[number_of_layer];

   ele_index_of_layer[0] = ele_index^1;
   for (size_t i = 1; i < number_of_layer - 1; i ++){
      ele_index_of_layer[i] = ele_index_of_layer[i - 1] >> 1;
      ele_index_of_layer[i] ^= 1;
   }
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      size_t size;
      if (i == 0){size = size_per_item;} else {size = SHA256_DIGEST_LENGTH;}
      gen_params(item_of_layer[i], size, N, logt, d, params, pir_params);
      
      alloc_server.construct(server + i, params, pir_params);
      alloc_client.construct(client + i, params, pir_params);
      GaloisKeys galois_keys = client[i].generate_galois_keys();
      server[i].set_galois_key(0, galois_keys);
      server[i].set_database(move(db[i]), item_of_layer[i], size);
      server[i].preprocess_database();
   }
   
   auto time_init_e = high_resolution_clock::now();
   auto time_init_us = duration_cast<microseconds>(time_init_e - time_init_s).count();
   cout << "____________________\n";
   std::cout << "Initialization total take time " << time_init_us / 1000 << " ms\n";
   cout << "____________________\n";

// Query
   auto time_query_s = high_resolution_clock::now();
   
   uint64_t index = client[0].get_fv_index(ele_index, size_per_item);
   uint64_t offset = client[0].get_fv_offset(ele_index, size_per_item); 
   PirQuery query = client[0].generate_query(index);
   PirQuery query_proof[number_of_layer - 1];
   for (size_t i = 0; i < number_of_layer - 1; i ++)
   {
      size_t size;
      if (i == 0){size = size_per_item;} else {size = SHA256_DIGEST_LENGTH;} 
      uint64_t index_proof = client[i].get_fv_index(ele_index_of_layer[i], size);   

      query_proof[i] = client[i].generate_query(index_proof);
   }

   auto time_query_e = high_resolution_clock::now();
   auto time_query_us = duration_cast<microseconds>(time_query_e - time_query_s).count();
   cout << "____________________\n";
   std::cout << "Query total take time " << time_query_us / 1000 << " ms\n";
   cout << "____________________\n";

// Reply
   auto time_reply_s = high_resolution_clock::now();

   PirReply reply = server[0].generate_reply(query, 0);
   PirReply reply_proof[number_of_layer - 1];
   for (size_t i = 0; i < number_of_layer - 1; i ++)
   {
      reply_proof[i] = server[i].generate_reply(query_proof[i], 0);
   }

   auto time_reply_e = high_resolution_clock::now();
   auto time_reply_us = duration_cast<microseconds>(time_reply_e - time_reply_s).count();
   cout << "____________________\n";
   std::cout << "Reply total take time " << time_reply_us / 1000 << " ms\n";
   cout << "____________________\n";

//Decode
   auto time_decode_s = high_resolution_clock::now();
   unsigned char * client_query_proof[number_of_layer - 1];
   unsigned char  client_query_db[size_per_item];
   uint8_t elems[N *logt / 8];
   Plaintext result = client[0].decode_reply(reply);
   coeffs_to_bytes(logt, result, elems, (N * logt) / 8);
   for (size_t i = 0; i < size_per_item; i ++){
      client_query_db[i] = elems[size_per_item * offset + i];
   }


   for (size_t i = 0; i < number_of_layer - 1; i ++){
      size_t size;
      if (i == 0){size = size_per_item;} else {size = SHA256_DIGEST_LENGTH;} 
// cout << "i = " << i << "\n";
// cout << "size = " << size << "\n";
// cout << "query index = " << ele_index_of_layer[i] << "\n";
      uint64_t index_proof = client[i].get_fv_index(ele_index_of_layer[i], size);   
      uint64_t offset_proof = client[i].get_fv_offset(ele_index_of_layer[i], size); 

      Plaintext result_proof = client[i].decode_reply(reply_proof[i]);
      uint8_t elems_proof[N * logt / 8];
      coeffs_to_bytes(logt, result_proof, elems_proof, (N * logt) / 8);

      bool flag = true;
      for (size_t j = 0 ; j < size; j ++){
         if (elems_proof[(offset_proof * size) + j] != db_copy[i].get()[(ele_index_of_layer[i] * size) + j]) {
            cout << j << "wrong\n";
            flag = false;
         }
      }
      // if (flag){
      //    cout << "layer " << i << " query correct.\n";
      // }
      
      client_query_proof[i] = (unsigned char *)malloc(sizeof(unsigned char) * size);
      for (size_t j = 0; j < size; j ++){
         client_query_proof[i][j] = elems_proof[offset_proof * size + j];
      }
   }
   
   auto time_decode_e = high_resolution_clock::now();
   auto time_decode_us = duration_cast<microseconds>(time_decode_e - time_decode_s).count();
   cout << "____________________\n";
   std::cout << "Decode total take time " << time_decode_us / 1000 << " ms\n";
   cout << "____________________\n";

// // Verify

   auto time_verify_s = high_resolution_clock::now();
   unsigned char data[2 * size_per_item];
   unsigned char hash[SHA256_DIGEST_LENGTH];
   for (size_t i = 0; i < size_per_item; i ++){
      if (ele_index % 2 == 0){
         data[i] = client_query_db[i];
         data[i + size_per_item] = client_query_proof[0][i];
      }
      else{
         data[i] = client_query_proof[0][i];
         data[i + size_per_item] = client_query_db[i];
      }
   }
   sha256(data, 2 * size_per_item, hash);
   bool flag = true;
   for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i ++){
      if (db_copy[1][ele_index_of_layer[1] * SHA256_DIGEST_LENGTH + i] != hash[i]){
         flag = false;
      }
   }
   // if (flag){
   //    cout << "layer 1 hash correct\n";
   // }
   // else{
   //    cout << "layer 1 hash wrong\n";
   //    cout << "should be hash(";
   //    unsigned char data2[2 * size_per_item];
   //    unsigned char hash2[SHA256_DIGEST_LENGTH];
      
   //    if (ele_index % 2 == 0){
   //       for (size_t i = 0; i < size_per_item; i++){
   //          cout << (int)db_copy[0][ele_index * size_per_item + i] << " ";
   //          data2[i] = db_copy[0][ele_index * size_per_item + i];
   //       }
   //       cout << "|";
   //       for (size_t i = 0; i < size_per_item; i ++){
   //          cout << (int)db_copy[0][ele_index_of_layer[0] * size_per_item + i] << " ";
   //          data2[i + size_per_item] = db_copy[0][ele_index_of_layer[0] * size_per_item + i];
   //       }
   //    }
   //    else{
   //       for (size_t i = 0; i < size_per_item; i ++){
   //          cout << (int)db_copy[0][ele_index_of_layer[0] * size_per_item + i] << " ";
   //          data2[i] = db_copy[0][ele_index_of_layer[0] * size_per_item + i];
   //       }
   //       cout << "|";
   //       for (size_t i = 0; i < size_per_item; i++){
   //          cout << (int)db_copy[0][ele_index * size_per_item + i] <<  " ";
   //          data2[i + size_per_item] = db_copy[0][ele_index * size_per_item + i];
   //       }
   //    }
   //    cout << ")\n =";
   //    for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i ++){
   //       cout << (int)db_copy[1][ele_index_of_layer[1] * SHA256_DIGEST_LENGTH + i] << " ";
   //    }
   //    sha256(data2, 2 * size_per_item, hash2);
      // cout << "\n = ";
      // for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i ++){
      //    cout << (int) hash2[i] << " ";
      // }

      // cout << "\n actually is hash(";
      // for (size_t i = 0; i < 2 * size_per_item; i ++){
      //    cout << (int)data[i] << " ";
      //    if (i == size_per_item - 1){cout << "|";}
      // }
      // cout << ")\n =";
      // for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i ++){
      //    cout << (int)hash[i] << " ";
      // }
   // }


   unsigned char data2[2 * SHA256_DIGEST_LENGTH];
   for (size_t i = 1; i < number_of_layer - 1; i++){
      for (size_t j = 0; j < SHA256_DIGEST_LENGTH; j ++){
         if (ele_index_of_layer[i] % 2 == 1){
            data2[j] = hash[j];
            data2[j + SHA256_DIGEST_LENGTH] = client_query_proof[i][j];
         }
         else{
            data2[j] = client_query_proof[i][j];
            data2[j + SHA256_DIGEST_LENGTH] = hash[j];
         }
         sha256(data2, 2 * SHA256_DIGEST_LENGTH, hash);
         // for (size_t j = 0; j < SHA256_DIGEST_LENGTH; j ++){
         //    if (db_copy[i + 1][ele_index_of_layer[i + 1] * SHA256_DIGEST_LENGTH + j] != hash[j]){
         //       flag = false;
         //       cout << "layer " << i+1 << " hash wrong\n";
         //    }
         // }
      }
   }
   
   auto time_verify_e = high_resolution_clock::now();
   auto time_verify_us = duration_cast<microseconds>(time_verify_e - time_verify_s).count();
   
   std::cout << "Verify total take time " << time_verify_us / 1000 << " ms\n";
   cout << "____________________\n";
   

   // for (size_t i = 0; i < SHA256_DIGEST_LENGTH; i ++){
   //    if (db[number_of_layer - 1][i] != hash[i]){
   //       cout << "wrong\n";
   //    }
   // }
   cout << "correct\n";

   size_t count = 0; 
   for (size_t i = 0; i < query.size(); i ++){
      count += query[i].size();
   }
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      for (size_t j = 0; j < query_proof[i].size(); j ++){
         count += query_proof[i][j].size();
      }
   }
   cout << "query size " << count * N * logt / 8 << " Bytes\n";
   
   size_t count2 = 0;
   count2 += reply.size();
   for (size_t i = 0; i < number_of_layer - 1; i ++){
      count2 += reply_proof[i].size();
   }
   cout << "reply size " << count2 * N * logt / 8 << " Bytes\n"; 

   

   return 0;
}