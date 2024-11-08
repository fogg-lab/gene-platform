// Import workers
import PyWorker from './PyWorker.worker';
import RWorker from './RWorker.worker';
import RustWorker from './RustWorker.worker';

class WorkerManager {
    constructor() {
        this.workers = {
            py: new PyWorker(),
            r: new RWorker(),
            rust: new RustWorker()
        };

        // Initialize error handling for each worker
        Object.values(this.workers).forEach(worker => {
            worker.onerror = (error) => {
                console.error('Worker error:', error);
            };
        });
    }

    async runTask(workerType, action, data) {
        return new Promise((resolve, reject) => {
            const worker = this.workers[workerType];
            
            if (!worker) {
                reject(new Error(`Invalid worker type: ${workerType}`));
                return;
            }

            const messageHandler = (event) => {
                if (event.data.status === 'success') {
                    worker.removeEventListener('message', messageHandler);
                    resolve(event.data.result);
                } else if (event.data.status === 'error') {
                    worker.removeEventListener('message', messageHandler);
                    reject(event.data.error);
                }
            };

            worker.addEventListener('message', messageHandler);

            try {
                worker.postMessage({ action, data });
            } catch (error) {
                worker.removeEventListener('message', messageHandler);
                reject(error);
            }
        });
    }

    terminateWorkers() {
        Object.values(this.workers).forEach(worker => {
            try {
                worker.terminate();
            } catch (error) {
                console.error('Error terminating worker:', error);
            }
        });
    }
}

const workerManager = new WorkerManager();
export default workerManager;
