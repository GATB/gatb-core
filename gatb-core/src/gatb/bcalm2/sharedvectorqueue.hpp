#pragma once
// from http://stackoverflow.com/questions/36762248/why-is-stdqueue-not-thread-safe

#include <queue>
#include <mutex>
#include <condition_variable>

    template <typename T>
    class SharedVectorQueue
    {
    public:
        SharedVectorQueue(int nb_elts) : nb_elts(nb_elts) {
        
            std::vector<std::mutex> list(nb_elts);
            mutex_.swap(list);

            queue_.resize(nb_elts);

            std::vector<std::condition_variable> cond_2(nb_elts);
            cond_.swap(cond_2);
        };

        ~SharedVectorQueue(){
            for (int i = 0; i < nb_elts; i++)
            {
                //std::deque<T>().swap(queue_[i]); 
                std::vector<T>().swap(queue_[i]); 
            }
        };
        SharedVectorQueue(const SharedVectorQueue&) = delete;

        bool pop_immediately(int i, T& item)
        {
            if (queue_[i].empty())
            {
                queue_[i].shrink_to_fit();
                return false;
            }
             // std::vector queue
            item = queue_[i].back();
            queue_[i].pop_back();
            
            /*
            item = queue_[i].front();
            queue_[i].pop_front();
            */
            return true;
        }

        void push_back(int i, const T& item)
        {
            std::unique_lock<std::mutex> mlock(mutex_[i]);
            queue_[i].push_back(item);
            mlock.unlock();     // unlock before notificiation to minimize mutex con
            cond_[i].notify_one(); // notify one waiting thread

        }
        void push_back(int i, T&& item)
        {
            std::unique_lock<std::mutex> mlock(mutex_[i]);
            queue_[i].push_back(std::move(item));
            mlock.unlock();     // unlock before notificiation to minimize mutex con
            cond_[i].notify_one(); // notify one waiting thread

        }


        int size(int i)
        {
            std::unique_lock<std::mutex> mlock(mutex_[i]);
            int size = queue_[i].size();
            mlock.unlock();
            return size;
        }    
        
        bool empty(int i)
        {
            return queue_[i].empty();
        }


    private:
        int nb_elts;
        std::vector<std::vector<T>> queue_;
        std::vector<std::mutex> mutex_;
        std::vector<std::condition_variable> cond_;
    }; 


         

    

    
