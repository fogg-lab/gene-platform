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
    }

    async runTask(workerType, action, data) {
      return new Promise((resolve, reject) => {
        const worker = this.workers[workerType];
        worker.onmessage = (event) => {
          if (event.data.status === 'success') {
            resolve(event.data.result);
          } else {
            reject(event.data.error);
          }
        };
        worker.postMessage({ action, data });
      });
    }

    terminateWorkers() {
      Object.values(this.workers).forEach(worker => worker.terminate());
    }
  }

  export default new WorkerManager();
