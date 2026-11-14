// TwistNTurn Service Worker
const CACHE_NAME = 'twistnturn-v3';

// Install event: prepare cache
self.addEventListener('install', event => {
    self.skipWaiting();
});

// Activate event: clean up outdated caches and take control
self.addEventListener('activate', event => {
    event.waitUntil(
        caches.keys().then(keys => Promise.all(
            keys.filter(k => k !== CACHE_NAME).map(k => caches.delete(k))
        )).then(() => self.clients.claim())
    );
});

// Fetch event:
// - Network-first for navigation and index.html (guarantees fast updates)
// - Cache-first with network fallback for fingerprinted assets
self.addEventListener('fetch', event => {
    if (event.request.method !== 'GET') return;

    const url = new URL(event.request.url);

    // Only cache requests from our origin
    if (url.origin !== self.location.origin) return;

    // Navigation / HTML: Network-first, fallback to cache
    if (event.request.mode === 'navigate' || url.pathname.endsWith('/') || url.pathname.endsWith('.html')) {
        event.respondWith(
            fetch(event.request)
                .then(networkResponse => {
                    if (networkResponse && networkResponse.status === 200) {
                        const copy = networkResponse.clone();
                        caches.open(CACHE_NAME).then(cache => cache.put(event.request, copy));
                    }
                    return networkResponse;
                })
                .catch(() => caches.match(event.request).then(cached => cached || caches.match('./')))
        );
        return;
    }

    // Static assets (_framework, css, images): Cache-first
    event.respondWith(
        caches.match(event.request).then(cachedResponse => {
            if (cachedResponse) {
                return cachedResponse;
            }
            return fetch(event.request).then(networkResponse => {
                if (networkResponse && networkResponse.status === 200) {
                    const copy = networkResponse.clone();
                    caches.open(CACHE_NAME).then(cache => cache.put(event.request, copy));
                }
                return networkResponse;
            });
        })
    );
});

// Message listener for skip waiting or cache busting
self.addEventListener('message', event => {
    if (event.data) {
        if (event.data.type === 'SKIP_WAITING') {
            self.skipWaiting();
        } else if (event.data.type === 'CLEAR_CACHE') {
            event.waitUntil(
                caches.keys().then(keys => Promise.all(keys.map(k => caches.delete(k))))
            );
        }
    }
});

